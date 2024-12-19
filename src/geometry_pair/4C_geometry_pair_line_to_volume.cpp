// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_volume.hpp"

#include "4C_geometry_pair_constants.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
GEOMETRYPAIR::GeometryPairLineToVolume<ScalarType, Line, Volume>::GeometryPairLineToVolume(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const std::shared_ptr<GEOMETRYPAIR::LineTo3DEvaluationData>& line_to_3d_evaluation_data)
    : GeometryPair(element1, element2), line_to_3d_evaluation_data_(line_to_3d_evaluation_data)
{
  // For the current implementation, the line element has to be on the same processor as the pair
  // object. This is because the tracking vector in LineTo3DEvaluationData is only local and we
  // need this vector for segmentation e.t.c.
  int myrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (element1->owner() != myrank)
    FOUR_C_THROW(
        "The GeometryPairLineToVolume pair has to be on the same processor as the line element! "
        "Currently the pair is on rank %d, the line element on %d!",
        myrank, element1->owner());
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<ScalarType, Line, Volume>::project_point_to_other(
    const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
    const ElementData<Volume, ScalarType>& element_data_volume,
    Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result) const
{
  GEOMETRYPAIR::project_point_to_volume<ScalarType, Volume>(
      point, element_data_volume, xi, projection_result);
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<ScalarType, Line, Volume>::intersect_line_with_surface(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Volume, ScalarType>& element_data_volume, const unsigned int& fixed_parameter,
    const double& fixed_value, ScalarType& eta, Core::LinAlg::Matrix<3, 1, ScalarType>& xi,
    ProjectionResult& projection_result) const
{
  // Check the input parameters.
  {
    if (Volume::geometry_type_ == DiscretizationTypeGeometry::hexahedron && fixed_parameter > 2)
      FOUR_C_THROW(
          "Fixed_parameter in IntersectLineWithVolume has to be smaller than 3 with a hexahedron "
          "element.");
    else if (Volume::element_dim_ != 3)
      FOUR_C_THROW("Wrong DiscretizationTypeGeometry type given.");
    else if (fixed_parameter > 3)
      FOUR_C_THROW("fixed_parameter in IntersectLineWithVolume can be 3 at maximum.");
  }

  // Initialize data structures
  // Point on line.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_line;
  Core::LinAlg::Matrix<3, 1, ScalarType> dr_line;

  // Point on volume.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_volume;
  Core::LinAlg::Matrix<3, 3, ScalarType> dr_volume;

  // Residuum.
  Core::LinAlg::Matrix<4, 1, ScalarType> residuum;
  Core::LinAlg::Matrix<4, 1, ScalarType> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);

  // Jacobian / inverse.
  Core::LinAlg::Matrix<4, 4, ScalarType> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  {
    // Local Newton iteration.
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Get the point coordinates on the line and volume.
      evaluate_position<Line>(eta, element_data_line, r_line);
      evaluate_position<Volume>(xi, element_data_volume, r_volume);

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$ and $xi(i) - value = R_{surf}$
      J_J_inv.put_scalar(0.);
      residuum.put_scalar(0.);
      for (unsigned int i = 0; i < 3; i++)
      {
        residuum(i) = r_volume(i) - r_line(i);
      }
      if (fixed_parameter < 3)
      {
        residuum(3) = xi(fixed_parameter) - fixed_value;
        J_J_inv(3, fixed_parameter) = 1.;
      }
      else
      {
        for (unsigned int i = 0; i < 3; i++)
        {
          residuum(3) += xi(i);
          J_J_inv(3, i) = 1.;
        }
        residuum(3) -= fixed_value;
      }

      // Check if tolerance is fulfilled.
      if (residuum.norm2() < Constants::local_newton_res_tol &&
          delta_xi.norm2() < Constants::projection_xi_eta_tol)
      {
        // Check if the parameter coordinates are valid.
        if (valid_parameter_1d(eta) && valid_parameter_3d<Volume>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.norm2() > Constants::local_newton_res_max) break;

      // Get the positional derivatives.
      evaluate_position_derivative1<Line>(eta, element_data_line, dr_line);
      evaluate_position_derivative1<Volume>(xi, element_data_volume, dr_volume);

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J_J_inv(i, j) = dr_volume(i, j);
        }
        J_J_inv(i, 3) = -dr_line(i);
      }

      // Solve the linearized system.
      if (Core::LinAlg::solve_linear_system_do_not_throw_error_on_zero_determinant_scaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        eta -= delta_xi(3);
        for (unsigned int i = 0; i < 3; i++) xi(i) -= delta_xi(i);

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<ScalarType, Line, Volume>::intersect_line_with_other(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Volume, ScalarType>& element_data_volume,
    std::vector<ProjectionPoint1DTo3D<ScalarType>>& intersection_points,
    const ScalarType& eta_start, const Core::LinAlg::Matrix<3, 1, ScalarType>& xi_start) const
{
  // Get number of faces for this volume and create a vector with the indices of the faces, so all
  // surfaces of the volume can be checked for an intersection with the line.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  if (Volume::geometry_type_ == DiscretizationTypeGeometry::hexahedron)
  {
    n_faces = 6;
    face_fixed_parameters = {0, 0, 1, 1, 2, 2};
    face_fixed_values = {-1., 1., -1., 1., -1., 1.};
  }
  else if (Volume::geometry_type_ == DiscretizationTypeGeometry::tetraeder)
  {
    n_faces = 4;
    face_fixed_parameters = {0, 1, 2, 3};
    face_fixed_values = {0., 0., 0., 1.};
  }
  else
  {
    FOUR_C_THROW("Wrong DiscretizationTypeGeometry given!");
  }

  // Clear the input vector.
  intersection_points.clear();
  intersection_points.reserve(n_faces);

  // Create variables.
  ScalarType eta;
  Core::LinAlg::Matrix<3, 1, ScalarType> xi;
  ProjectionResult intersection_found;

  // Try to intersect the beam with each face.
  for (unsigned int i = 0; i < n_faces; i++)
  {
    // Set starting values.
    xi = xi_start;
    eta = eta_start;

    // Intersect the line with the surface.
    intersect_line_with_surface(element_data_line, element_data_volume, face_fixed_parameters[i],
        face_fixed_values[i], eta, xi, intersection_found);

    // If a valid intersection is found, add it to the output vector.
    if (intersection_found == ProjectionResult::projection_found_valid)
    {
      intersection_points.push_back(ProjectionPoint1DTo3D<ScalarType>(eta, xi));
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Volume>
void GEOMETRYPAIR::project_point_to_volume(const Core::LinAlg::Matrix<3, 1, ScalarType>& point,
    const ElementData<Volume, ScalarType>& element_data_volume,
    Core::LinAlg::Matrix<3, 1, ScalarType>& xi, ProjectionResult& projection_result)
{
  // Point on volume.
  Core::LinAlg::Matrix<3, 1, ScalarType> r_volume;

  // Jacobian / inverse.
  Core::LinAlg::Matrix<3, 3, ScalarType> J_J_inv;

  // Increment of xi.
  Core::LinAlg::Matrix<3, 1, ScalarType> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);

  // Residuum.
  Core::LinAlg::Matrix<3, 1, ScalarType> residuum;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Get the point coordinates on the volume.
      GEOMETRYPAIR::evaluate_position<Volume>(xi, element_data_volume, r_volume);

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$
      residuum = r_volume;
      residuum -= point;

      // Check if tolerance is fulfilled.
      if (residuum.norm2() < Constants::local_newton_res_tol &&
          delta_xi.norm2() < Constants::projection_xi_eta_tol)
      {
        if (valid_parameter_3d<Volume>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.norm2() > Constants::local_newton_res_max) break;

      // Get the jacobian.
      GEOMETRYPAIR::evaluate_position_derivative1<Volume>(xi, element_data_volume, J_J_inv);

      // Solve the linearized system.
      if (Core::LinAlg::solve_linear_system_do_not_throw_error_on_zero_determinant_scaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        xi -= delta_xi;

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_wedge6>;

FOUR_C_NAMESPACE_CLOSE
