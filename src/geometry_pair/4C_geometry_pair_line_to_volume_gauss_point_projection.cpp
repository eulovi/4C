// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometry_pair_line_to_volume_gauss_point_projection.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_line_projection.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_geometry_pair_utility_classes.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line,
    Volume>::GeometryPairLineToVolumeGaussPointProjection(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const std::shared_ptr<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
    : GeometryPairLineToVolume<ScalarType, Line, Volume>(element1, element2, evaluation_data)
{
  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->element1()->id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->get_gauss_point_projection_tracker();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points = this->line_to_3d_evaluation_data_->get_number_of_gauss_points();
    std::vector<bool> new_tracking_vector;
    new_tracking_vector.resize(n_gauss_points, false);
    projection_tracker[line_element_id] = new_tracking_vector;
  }
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line,
    Volume>::pre_evaluate(const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Volume, ScalarType>& element_data_volume,
    std::vector<LineSegment<ScalarType>>& segments) const
{
  // Call the pre_evaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line, Volume>>::pre_evaluate(this,
      element_data_line, element_data_volume, segments);
}


/**
 *
 */
template <typename ScalarType, typename Line, typename Volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line, Volume>::evaluate(
    const ElementData<Line, ScalarType>& element_data_line,
    const ElementData<Volume, ScalarType>& element_data_volume,
    std::vector<LineSegment<ScalarType>>& segments) const
{
  // Call the pre_evaluate method of the general Gauss point projection class.
  LineTo3DGaussPointProjection<
      GeometryPairLineToVolumeGaussPointProjection<ScalarType, Line, Volume>>::evaluate(this,
      element_data_line, element_data_volume, segments);
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_wedge6>;

FOUR_C_NAMESPACE_CLOSE
