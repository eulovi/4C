// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_INTERACTION_PAIR_NITSCHE_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_INTERACTION_PAIR_NITSCHE_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_embeddedmesh_interaction_pair.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Constraints::EmbeddedMesh
{
  class SolidToSolidNitscheManager;

  template <typename Interface, typename ParentInterface, typename Background>
  class SurfaceToBackgroundCouplingPairNitsche : public SolidInteractionPair
  {
   public:
    /**
     * \brief Standard Constructor
     */
    SurfaceToBackgroundCouplingPairNitsche(std::shared_ptr<Core::Elements::Element> element1,
        Core::Elements::Element* element2,
        Constraints::EmbeddedMesh::EmbeddedMeshParams& params_ptr,
        std::shared_ptr<Cut::CutWizard>& cutwizard_ptr,
        std::vector<std::shared_ptr<Cut::BoundaryCell>>& boundary_cells);

    /**
     * \brief Destructor.
     */
    ~SurfaceToBackgroundCouplingPairNitsche() override {};

    void get_pair_visualization(Core::IO::VisualizationData& visualization_data,
        std::shared_ptr<Core::LinAlg::Vector<double>> vector_to_visualize,
        const Constraints::EmbeddedMesh::SolidToSolidCouplingManager* coupling_manager,
        std::shared_ptr<std::unordered_set<int>> interface_tracker) override
    {
      FOUR_C_THROW("This function is not implemented for a Nitsche pair");
    }

    //! @name Evaluation methods
    /**
     * \brief Evaluate the global matrices and vectors resulting from the coupling method.
     * @param discret (in) Discretization, used to get the boundary layer GIDs.
     * @param coupling_manager (in) Manager for the specific type of coupling technique
     * @param global_matrices (in) A map of sparse matrices resulting from the coupling method.
     * @param global_vectors (in) A map of FEVectors resulting from the coupling method.
     */
    void evaluate_and_assemble_coupling_contributions(const Core::FE::Discretization& discret,
        const Constraints::EmbeddedMesh::SolidToSolidCouplingManager* coupling_manager,
        std::map<std::string, std::shared_ptr<Core::LinAlg::SparseMatrix>>& global_matrices,
        std::map<std::string, std::shared_ptr<Core::LinAlg::FEVector<double>>>& global_vectors)
        override;

    /**
     * \brief Evaluate the global matrices and vectors resulting from Nitsche coupling.
     * @param discret (in) Discretization, used to get the interface GIDs.
     * @param nitsche_manager (in) Nitsche manager, used to get the Lagrange multiplier GIDs.
     * TODO: add description of variables
     */
    void evaluate_and_assemble_nitsche_contributions(const Core::FE::Discretization& discret,
        const Constraints::EmbeddedMesh::SolidToSolidNitscheManager* nitsche_manager,
        Core::LinAlg::SparseMatrix& global_penalty_interface,
        Core::LinAlg::SparseMatrix& global_penalty_background,
        Core::LinAlg::SparseMatrix& global_penalty_interface_background,
        Core::LinAlg::SparseMatrix& global_nitsche_parent,
        Core::LinAlg::SparseMatrix& global_nitsche_background,
        Core::LinAlg::SparseMatrix& global_nitsche_parent_background,
        Core::LinAlg::FEVector<double>& global_penalty_constraint,
        Core::LinAlg::FEVector<double>& global_nitsche_constraint,
        double& nitsche_average_weight_param);

    /**
     * \brief Set the Gauss rule over the interface for element1_ and element2_.
     */
    void set_gauss_rule_for_interface_and_background();

    /**
     * \brief Update the current displacement of the interface and background elements
     */
    void set_current_element_position(Core::FE::Discretization const& discret,
        const Core::LinAlg::Vector<double>& displacement_vector) override;

   private:
    /**
     * \brief Evaluate the local penalty contributions of the Nitsche method for this coupling
     * element pair.
     */
    void evaluate_penalty_contributions_nitsche(
        Core::LinAlg::Matrix<Interface::n_dof_, Interface::n_dof_, double>&
            local_stiffness_penalty_interface,
        Core::LinAlg::Matrix<Background::n_dof_, Background::n_dof_, double>&
            local_stiffness_penalty_background,
        Core::LinAlg::Matrix<Interface::n_dof_, Background::n_dof_, double>&
            local_stiffness_penalty_interface_background,
        Core::LinAlg::Matrix<Interface::n_dof_ + Background::n_dof_, 1, double>& local_constraint);

    void evaluate_stress_contributions_nitsche(const Core::FE::Discretization& discret,
        Core::LinAlg::Matrix<ParentInterface::n_dof_, ParentInterface::n_dof_, double>&
            local_stiffness_nitsche_parent,
        Core::LinAlg::Matrix<Background::n_dof_, Background::n_dof_, double>&
            local_stiffness_nitsche_background,
        Core::LinAlg::Matrix<ParentInterface::n_dof_, Background::n_dof_, double>&
            local_stiffness_nitsche_parent_background,
        Core::LinAlg::Matrix<ParentInterface::n_dof_ + Background::n_dof_, 1, double>&
            local_constraint_nitsche,
        double& nitsche_average_weight_param);

    //! @name Visualization methods
    /*!
    \brief Get the Gauss points of element2_ after cut for visualization.
    */
    void get_projected_gauss_rule_in_cut_element(
        Core::IO::VisualizationData& cut_element_integration_points_visualization_data) override;

    /*!
    \brief Get the Gauss points of element1_ and element2_ for the evaluation of mortar matrices.
    */
    void get_projected_gauss_rule_on_interface(
        Core::IO::VisualizationData& background_integration_points_visualization_data,
        Core::IO::VisualizationData& interface_integration_points_visualization_data) override;

    //! Get parent element of first element
    Core::Elements::Element& parent_element_1() const { return *parent_element1_; };

    //! Kinematics of parent element
    FourC::Solid::KinemType parent_element1_kinem_;

    //! Kinematics of background element
    FourC::Solid::KinemType element2_kinem_;

    //! Pointer to the parent element of element_1()
    Core::Elements::Element* parent_element1_ = nullptr;

    //! Initial nodal positions (and tangents) of the interface element.
    GeometryPair::ElementData<Interface, double> ele1pos_;

    //! Initial nodal positions (and tangents) of the parent element of the interface element.
    GeometryPair::ElementData<ParentInterface, double> ele1parentpos_;

    //! Initial nodal positions (and tangents) of the background element.
    GeometryPair::ElementData<Background, double> ele2pos_;

    //! Displacements of the interface element.
    GeometryPair::ElementData<Interface, double> ele1dis_;

    //! Displacements of the parent of the interface element.
    GeometryPair::ElementData<ParentInterface, double> ele1parentdis_;

    //! Displacements of the background element.
    GeometryPair::ElementData<Background, double> ele2dis_;

    //! integration rule over the interface for element1_ and element2_
    std::vector<std::tuple<Core::LinAlg::Matrix<2, 1>, Core::LinAlg::Matrix<3, 1>, double>>
        interface_integration_points_;
  };
}  // namespace Constraints::EmbeddedMesh

FOUR_C_NAMESPACE_CLOSE

#endif
