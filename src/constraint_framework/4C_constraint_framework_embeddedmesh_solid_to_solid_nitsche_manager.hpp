// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_SOLID_TO_SOLID_NITSCHE_MANAGER_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_SOLID_TO_SOLID_NITSCHE_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_constraint_framework_embeddedmesh_params.hpp"
#include "4C_constraint_framework_embeddedmesh_solid_to_solid_coupling_manager.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_fevector.hpp"
#include "4C_utils_exceptions.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
class Map;

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Core::IO
{
  class VisualizationManager;
}

namespace Solid
{
  namespace TimeInt
  {
    class BaseDataGlobalState;
  }
}  // namespace Solid

namespace LinAlg
{
  class SparseMatrix;
}  // namespace LinAlg

namespace Constraints::EmbeddedMesh
{
  class SolidInteractionPair;

  class SolidToSolidNitscheManager : public SolidToSolidCouplingManager
  {
   public:
    /**
     * \brief Standard Constructor
     *
     * @param discret (in) Pointer to the discretization.
     * @param displacement_vector (in) global displacement vector.
     * @param embedded_mesh_coupling_params (in) embedded mesh coupling parameters
     * @param visualization_manager (in) visualization manager
     */
    SolidToSolidNitscheManager(std::shared_ptr<Core::FE::Discretization>& discret,
        const Core::LinAlg::Vector<double>& displacement_vector,
        Constraints::EmbeddedMesh::EmbeddedMeshParams& embedded_mesh_coupling_params,
        std::shared_ptr<Core::IO::VisualizationManager> visualization_manager);

    /**
     * \brief This method builds the global maps for the global node
     * @param displacement_vector (in) global displacement vector.
     */
    void setup(const Core::LinAlg::Vector<double>& displacement_vector) override;

    /**
     * \brief Evaluate Nitsche coupling contributions on all pairs and assemble them into the
     * global matrices.
     * @param displacement_vector (in) global displacement vector.
     */
    void evaluate_global_coupling_contributions(
        const Core::LinAlg::Vector<double>& displacement_vector) override;

    /**
     *
     */
    void add_global_force_stiffness_contributions(Solid::TimeInt::BaseDataGlobalState& data_state,
        std::shared_ptr<Core::LinAlg::SparseMatrix> stiff,
        std::shared_ptr<Core::LinAlg::Vector<double>> force) const override;

    /**
     * \brief Scale contributions from the Nitsche method with the penalty parameter and weighting
     * average parameter
     */
    void scale_penalty_contributions() const;

    /**
     * \brief Write output obtained in the embedded mesh
     */
    void write_output(double time, int timestep_number) override;

    /**
     * \brief Obtain the energy contribution of the embedded mesh method
     */
    double get_energy() const override
    {
      FOUR_C_THROW("get_energy() not implemented for the Nitsche strategy");
    }

   protected:
    /**
     * \brief Throw an error if setup was not called on the object prior to this function call.
     */
    inline void check_setup() const override
    {
      if (!is_setup_) FOUR_C_THROW("Setup not called on SolidToSolidNitscheManager!");
    }

    /**
     * \brief Throw an error if the global maps were not build.
     */
    inline void check_global_maps() const override
    {
      if (!is_global_maps_build_)
        FOUR_C_THROW("Global maps are not build in SolidToSolidNitscheManager!");
    }

    /**
     * \brief Check if this node is in a parent element
     */
    bool is_parent_node(Core::Nodes::Node const& node);

   private:
    /**
     * \brief Calculate the maps for the solid interface and background dofs. The calculated
     * maps are used to complete the Nitsche matrices.
     */
    void set_global_maps() override;

    std::vector<int> parent_node_ids_;

    //! Row map of the solid interface DOFs.
    std::shared_ptr<Core::LinAlg::Map> interface_dof_rowmap_;

    //! Row map of the solid parent element of the interface DOFs.
    std::shared_ptr<Core::LinAlg::Map> parent_dof_rowmap_;

    //! Row map of both interface and solid background DOFs.
    std::shared_ptr<Core::LinAlg::Map> interface_and_background_dof_rowmap_;

    std::shared_ptr<Core::LinAlg::Map> parent_and_background_dof_rowmap_;

    //! Global contributions of the penalty term associated with the interface DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_penalty_interface_ = nullptr;

    //! Global contributions of the penalty term associated with the background DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_penalty_background_ = nullptr;

    //! Global contributions of the penalty term associated with both interface and background DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_penalty_interface_background_ = nullptr;

    //! Global contributions of the Nitsche method associated with the parent of the interface DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_nitsche_parent_ = nullptr;

    //! Global contributions of the Nitsche method associated with the background DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_nitsche_background_ = nullptr;

    //! Global contributions of the Nitsche method associated with both the parent of the interface
    //! and background DOFs
    std::shared_ptr<Core::LinAlg::SparseMatrix> global_nitsche_parent_background_ = nullptr;

    //! Global penalty constraint vector.
    std::shared_ptr<Core::LinAlg::FEVector<double>> global_penalty_constraint_ = nullptr;

    //! Global Nitsche constraint vector.
    std::shared_ptr<Core::LinAlg::FEVector<double>> global_nitsche_constraint_ = nullptr;
  };

}  // namespace Constraints::EmbeddedMesh
FOUR_C_NAMESPACE_CLOSE
#endif
