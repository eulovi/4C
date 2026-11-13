// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_constraint_framework_embeddedmesh_solid_to_solid_nitsche_manager.hpp"

#include "4C_constraint_framework_embeddedmesh_interaction_pair.hpp"
#include "4C_constraint_framework_embeddedmesh_params.hpp"
#include "4C_constraint_framework_embeddedmesh_solid_to_solid_utils.hpp"
#include "4C_cut_cutwizard.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_linalg_fevector.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
Constraints::EmbeddedMesh::SolidToSolidNitscheManager::SolidToSolidNitscheManager(
    std::shared_ptr<Core::FE::Discretization>& discret,
    const Core::LinAlg::Vector<double>& displacement_vector,
    Constraints::EmbeddedMesh::EmbeddedMeshParams& embedded_mesh_coupling_params,
    std::shared_ptr<Core::IO::VisualizationManager> visualization_manager)
    : SolidToSolidCouplingManager(discret, embedded_mesh_coupling_params, visualization_manager)
{
  // Get the indices of the DOFs related to the parent elements of the interface
  std::vector<int> parent_dofs_temp = Constraints::EmbeddedMesh::get_node_ids_parent_elements(
      info_background_interface_elements_, *discret_);
  parent_node_ids_.reserve(parent_dofs_temp.size());
  parent_node_ids_.insert(parent_node_ids_.end(), parent_dofs_temp.begin(), parent_dofs_temp.end());

  // Setup the solid to solid Nitsche manager
  Core::LinAlg::Vector<double> disp_col_vec(*discret_->dof_col_map());
  Core::LinAlg::export_to(displacement_vector, disp_col_vec);
  setup(disp_col_vec);
}

void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::setup(
    const Core::LinAlg::Vector<double>& displacement_vector)
{
  // Get the global ids of all mesh nodes on this rank
  std::vector<int> my_interface_nodes_gid;
  for (int i_node = 0; i_node < discret_->node_row_map()->num_my_elements(); i_node++)
  {
    Core::Nodes::Node const& node = *(discret_->l_row_node(i_node));
    if (Constraints::EmbeddedMesh::is_interface_node(*discret_, node))
      my_interface_nodes_gid.push_back(node.id());
  }

  // Create the maps for interface and background DOFs.
  set_global_maps();

  // Create the global coupling matrices related to the penalty contributions.
  global_penalty_interface_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *interface_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);
  global_penalty_background_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *background_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);
  global_penalty_interface_background_ =
      std::make_shared<Core::LinAlg::SparseMatrix>(*interface_and_background_dof_rowmap_, 100, true,
          true, Core::LinAlg::SparseMatrix::FE_MATRIX);

  // Create the global coupling matrices related to the Nitsche stress contributions.
  global_nitsche_parent_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *parent_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);
  global_nitsche_background_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *background_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);
  global_nitsche_parent_background_ = std::make_shared<Core::LinAlg::SparseMatrix>(
      *parent_and_background_dof_rowmap_, 100, true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);

  // Create the global constraint contributions to the force.
  global_penalty_constraint_ =
      std::make_shared<Core::LinAlg::FEVector<double>>(*interface_and_background_dof_rowmap_);
  global_nitsche_constraint_ =
      std::make_shared<Core::LinAlg::FEVector<double>>(*parent_and_background_dof_rowmap_);

  // Set flag for successful setup.
  is_setup_ = true;

  // Set current state
  set_state(displacement_vector);
}

void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::set_global_maps()
{
  // Get the dofs of the background and interface elements
  std::vector<int> interface_dofs;
  std::vector<int> parent_dofs;
  std::vector<int> background_dofs;
  std::vector<int> interface_and_background_dofs;
  std::vector<int> parent_and_background_dofs;
  for (int i_node = 0; i_node < discret_->node_row_map()->num_my_elements(); i_node++)
  {
    const Core::Nodes::Node* node = discret_->l_row_node(i_node);
    if (is_cut_node(*node))
      discret_->dof(node, background_dofs);
    else if (Constraints::EmbeddedMesh::is_interface_node(*discret_, *node))
      discret_->dof(node, interface_dofs);

    // std::cout << "set global maps" << std::endl;
    if (is_parent_node(*node))
    {
      // std::cout << "is parent node " << node->id() << std::endl;
      discret_->dof(node, parent_dofs);
      // for (int i = 0; i < parent_dofs.size(); i++)
      //   std::cout << parent_dofs[i] << ", ";
      // std::cout << std::endl;
    }
  }
  //
  // std::cout << "interface_dofs" << std::endl;
  // for (int i = 0; i < interface_dofs.size(); i++)
  //   std::cout << interface_dofs[i] << ", ";
  // std::cout << std::endl;
  //
  // std::cout << "parent_dofs" << std::endl;
  // for (int i = 0; i < parent_dofs.size(); i++)
  //   std::cout << parent_dofs[i] << ", ";
  // std::cout << std::endl;
  //
  // std::cout << "background_dofs" << std::endl;
  // for (int i = 0; i < background_dofs.size(); i++)
  //   std::cout << background_dofs[i] << ", ";
  // std::cout << std::endl;

  interface_and_background_dofs.reserve(interface_dofs.size() + background_dofs.size());
  interface_and_background_dofs.insert(
      interface_and_background_dofs.end(), interface_dofs.begin(), interface_dofs.end());
  interface_and_background_dofs.insert(
      interface_and_background_dofs.end(), background_dofs.begin(), background_dofs.end());

  std::sort(interface_and_background_dofs.begin(), interface_and_background_dofs.end());

  parent_and_background_dofs.reserve(parent_dofs.size() + background_dofs.size());
  parent_and_background_dofs.insert(
      parent_and_background_dofs.end(), parent_dofs.begin(), parent_dofs.end());
  parent_and_background_dofs.insert(
      parent_and_background_dofs.end(), background_dofs.begin(), background_dofs.end());

  std::sort(parent_and_background_dofs.begin(), parent_and_background_dofs.end());

  // std::cout << "parent dofs ids " << std::endl;
  // for (int i = 0; i < parent_dofs.size(); i++)
  //   std::cout << parent_dofs[i] << std::endl;

  // Create the interface and solid maps.
  interface_dof_rowmap_ = std::make_shared<Core::LinAlg::Map>(
      -1, interface_dofs.size(), interface_dofs.data(), 0, discret_->get_comm());
  parent_dof_rowmap_ = std::make_shared<Core::LinAlg::Map>(
      -1, parent_dofs.size(), parent_dofs.data(), 0, discret_->get_comm());
  background_dof_rowmap_ = std::make_shared<Core::LinAlg::Map>(
      -1, background_dofs.size(), background_dofs.data(), 0, discret_->get_comm());
  interface_and_background_dof_rowmap_ =
      std::make_shared<Core::LinAlg::Map>(-1, interface_and_background_dofs.size(),
          interface_and_background_dofs.data(), 0, discret_->get_comm());
  parent_and_background_dof_rowmap_ =
      std::make_shared<Core::LinAlg::Map>(-1, parent_and_background_dofs.size(),
          parent_and_background_dofs.data(), 0, discret_->get_comm());

  // Set flags for global maps.
  is_global_maps_build_ = true;
}


void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::evaluate_global_coupling_contributions(
    const Core::LinAlg::Vector<double>& displacement_vector)
{
  check_setup();
  check_global_maps();

  Core::LinAlg::Vector<double> disp_col_vec(*discret_->dof_col_map());
  Core::LinAlg::export_to(displacement_vector, disp_col_vec);
  set_state(disp_col_vec);

  // Reset the global data structures.
  global_penalty_interface_->put_scalar(0.);
  global_penalty_background_->put_scalar(0.);
  global_penalty_interface_background_->put_scalar(0.);

  global_nitsche_parent_->put_scalar(0.);
  global_nitsche_background_->put_scalar(0.);
  global_nitsche_parent_background_->put_scalar(0.);

  global_penalty_constraint_->put_scalar(0.);
  global_nitsche_constraint_->put_scalar(0.);

  std::map<std::string, std::shared_ptr<Core::LinAlg::SparseMatrix>> global_matrices = {
      {"global_penalty_interface", global_penalty_interface_},
      {"global_penalty_background", global_penalty_background_},
      {"global_penalty_interface_background", global_penalty_interface_background_},
      {"global_nitsche_parent", global_nitsche_parent_},
      {"global_nitsche_background", global_nitsche_background_},
      {"global_nitsche_parent_background", global_nitsche_parent_background_},
  };

  std::map<std::string, std::shared_ptr<Core::LinAlg::FEVector<double>>> global_vectors = {
      {"global_penalty_constraint", global_penalty_constraint_},
      {"global_nitsche_constraint", global_nitsche_constraint_},
  };

  for (auto& elepairptr : embedded_mesh_solid_pairs_)
  {
    elepairptr->evaluate_and_assemble_coupling_contributions(
        *discret_, this, global_matrices, global_vectors);
  }

  // Complete the global matrices.
  global_penalty_interface_->complete(*interface_dof_rowmap_, *interface_dof_rowmap_);
  global_penalty_background_->complete(*background_dof_rowmap_, *background_dof_rowmap_);
  global_penalty_interface_background_->complete(
      *interface_and_background_dof_rowmap_, *interface_and_background_dof_rowmap_);

  global_nitsche_parent_->complete(*parent_dof_rowmap_, *parent_dof_rowmap_);
  global_nitsche_background_->complete(*background_dof_rowmap_, *background_dof_rowmap_);
  global_nitsche_parent_background_->complete(
      *parent_and_background_dof_rowmap_, *parent_and_background_dof_rowmap_);

  // Complete the global vector.
  global_penalty_constraint_->complete();
  global_nitsche_constraint_->complete();

  // Scale the global penalty contributions before assembling into the final Nitsche formulation
  scale_penalty_contributions();
}

void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::scale_penalty_contributions() const
{
  // Get the penalty parameter.
  const double penalty_param = embedded_mesh_coupling_params_.nitsche_stabilization_param_;

  // Scale penalty contributions
  global_penalty_interface_->scale(penalty_param);
  global_penalty_background_->scale(penalty_param);
  global_penalty_interface_background_->scale(-penalty_param);

  global_penalty_constraint_->scale(penalty_param);
}

/**
 * here is the assemble of the contributions
 */
void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::
    add_global_force_stiffness_contributions(Solid::TimeInt::BaseDataGlobalState& data_state,
        std::shared_ptr<Core::LinAlg::SparseMatrix> stiff,
        std::shared_ptr<Core::LinAlg::Vector<double>> force) const
{
  check_setup();
  check_global_maps();

  if (stiff != nullptr)
  {
    stiff->add(*global_penalty_interface_, false, 1.0, 1.0);
    stiff->add(*global_penalty_background_, false, 1.0, 1.0);
    stiff->add(*global_penalty_interface_background_, false, 1.0, 1.0);
    stiff->add(*global_penalty_interface_background_, true, 1.0, 1.0);

    stiff->add(*global_nitsche_parent_, false, 1.0, 1.0);
    stiff->add(*global_nitsche_background_, false, 1.0, 1.0);
    stiff->add(*global_nitsche_parent_background_, false, 1.0, 1.0);
    stiff->add(*global_nitsche_parent_background_, true, 1.0, 1.0);
  }

  if (force != nullptr)
  {
    // Factor for right hand side (forces). 1 corresponds to forces being added to
    // the right hand side, -1 to the left hand side.
    const double rhs_factor = 1.0;

    // Collect force contributions in a global temp
    Core::LinAlg::Vector<double> global_penalty_temp(*discret_->dof_row_map());
    Core::LinAlg::export_to(*global_penalty_constraint_, global_penalty_temp);

    Core::LinAlg::Vector<double> global_nitsche_temp(*discret_->dof_row_map());
    Core::LinAlg::export_to(*global_nitsche_constraint_, global_nitsche_temp);

    // Add force contributions to global vector.
    force->update(rhs_factor, global_penalty_temp, 1.0);
    force->update(rhs_factor, global_nitsche_temp, 1.0);
  }
}


bool Constraints::EmbeddedMesh::SolidToSolidNitscheManager::is_parent_node(
    Core::Nodes::Node const& node)
{
  // Check if the id of this node is in the list of node ids.
  int node_id = node.id();
  return std::find(parent_node_ids_.begin(), parent_node_ids_.end(), node_id) !=
         parent_node_ids_.end();
}

void Constraints::EmbeddedMesh::SolidToSolidNitscheManager::write_output(
    double time, int timestep_number)
{
  // Clear the data of visualization manager
  visualization_manager_->clear_data();

  // Collect output of the physical positions of gauss points
  // where surface integrations will take place
  collect_output_integration_points();

  // Write output
  visualization_manager_->write_to_disk(time, timestep_number);
}

FOUR_C_NAMESPACE_CLOSE
