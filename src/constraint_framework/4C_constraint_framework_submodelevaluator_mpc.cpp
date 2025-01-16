// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_constraint_framework_submodelevaluator_mpc.hpp"

#include "4C_beam3_base.hpp"
#include "4C_constraint_framework_equation_mpc.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_mpc_rve.hpp"
#include "4C_io.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_OPEN
CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::RveMultiPointConstraintManager(
    std::shared_ptr<const Core::FE::Discretization> disc_ptr, Core::LinAlg::SparseMatrix* st_ptr)
{
  discret_ptr_ = disc_ptr;
  stiff_ptr_ = st_ptr;

  check_input();

  // Maps of relevant boundary nodesets of the rve
  std::map<std::string, const std::vector<int>*> rveBoundaryNodeIdMap;

  // Map of the corner node ids
  std::map<std::string, int> rveCornerNodeIdMap;

  //  Map the Node IDs to the respective rve boundary --> rveBoundaryNodeIdMap
  build_periodic_rve_boundary_node_map(rveBoundaryNodeIdMap);

  // Map the Node ids to the respective corner of the rve --> rveCornerNodeIdMap
  switch (rve_ref_type_)
  {
    case Inpar::RveMpc::RveReferenceDeformationDefinition::automatic:
    {
      build_periodic_rve_corner_node_map(rveBoundaryNodeIdMap, rveCornerNodeIdMap);
    }
    break;

    case Inpar::RveMpc::RveReferenceDeformationDefinition::manual:
    {
      if (rve_dim_ != Inpar::RveMpc::RveDimension::rve2d)
        FOUR_C_THROW("Manual Edge node definition is not implemented for 3D RVEs");

      // Read the reference points
      for (const auto& entry : point_periodic_rve_ref_conditions_)
      {
        const auto& str_id = entry->parameters().get<std::string>("POSITION");
        const auto* nodeInSet = entry->get_nodes();

        if (nodeInSet->size() > 1)
        {
          FOUR_C_THROW("There can only be a single node defined as a reference node");
        }
        rve_ref_node_map_[str_id.c_str()] = discret_ptr_->g_node(nodeInSet->data()[0]);

        Core::IO::cout(Core::IO::verbose)
            << "Map " << str_id.c_str() << " to node id " << (*nodeInSet)[0] << Core::IO::endl;
      }

      Core::IO::cout(Core::IO::verbose)
          << "Reference Points determined:"
          << "+--------------------------------------------------------------------+"
          << Core::IO::endl;
      for (const auto& elem : rve_ref_node_map_)
      {
        Core::IO::cout(Core::IO::verbose) << elem.first << ": " << elem.second->id() << ", ";
      }
      Core::IO::cout(Core::IO::verbose) << Core::IO::endl;

      // calculate the Reference vectors between Ref. poitns
      r_xmxp_[0] = rve_ref_node_map_["N2"]->x()[0] - rve_ref_node_map_["N1L"]->x()[0];
      r_xmxp_[1] = rve_ref_node_map_["N2"]->x()[1] - rve_ref_node_map_["N1L"]->x()[1];
      Core::IO::cout(Core::IO::verbose) << "RVE reference vector (X- ---> X+ ) : [" << r_xmxp_[0]
                                        << ", " << r_xmxp_[1] << "]" << Core::IO::endl;

      r_ymyp_[0] = rve_ref_node_map_["N4"]->x()[0] - rve_ref_node_map_["N1B"]->x()[0];
      r_ymyp_[1] = rve_ref_node_map_["N4"]->x()[1] - rve_ref_node_map_["N1B"]->x()[1];
      Core::IO::cout(Core::IO::verbose) << "RVE reference vector (Y- ---> Y+ ) : [" << r_ymyp_[0]
                                        << ", " << r_ymyp_[1] << "]" << Core::IO::endl;
    }
    break;
  }

  // Create a vector with all MPCs describing the periodic BCs
  if (surface_periodic_rve_conditions_.size() != 0 || line_periodic_rve_conditions_.size() != 0)
    build_periodic_mp_cs(rveBoundaryNodeIdMap, rveCornerNodeIdMap);


  // Add Linear Coupled Equation MPCs
  if (point_linear_coupled_equation_conditions_.size() != 0)
  {
    int nLinCe = build_linear_mp_cs();
    Core::IO::cout(Core::IO::verbose)
        << "Total number of linear coupled equations: " << nLinCe << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::check_input()
{
  if (Core::Communication::num_mpi_ranks(discret_ptr_->get_comm()) > 1)
    FOUR_C_THROW("periodic boundary conditions for RVEs are not implemented in parallel.");


  mpc_parameter_list_ = Global::Problem::instance()->rve_multi_point_constraint_params();

  rve_ref_type_ = Teuchos::getIntegralValue<Inpar::RveMpc::RveReferenceDeformationDefinition>(
      mpc_parameter_list_, "RVE_REFERENCE_POINTS");

  strategy_ = Teuchos::getIntegralValue<Inpar::RveMpc::EnforcementStrategy>(
      mpc_parameter_list_, "ENFORCEMENT");


  // Check the enforcement strategy
  switch (strategy_)
  {
    case Inpar::RveMpc::EnforcementStrategy::lagrangeMultiplier:
      FOUR_C_THROW("Constraint Enforcement via Lagrange Multiplier Methode is not impl.");

    case Inpar::RveMpc::EnforcementStrategy::penalty:
    {
      Core::IO::cout(Core::IO::minimal)
          << "Constraint enforcement strategy: Penalty method" << Core::IO::endl;
      const Teuchos::ParameterList& MpcList =
          Global::Problem::instance()->rve_multi_point_constraint_params();
      get_penalty_parameter_ptr() = MpcList.get<double>("PENALTY_PARAM");
      Core::IO::cout(Core::IO::verbose)
          << "Penalty weight used: " << get_penalty_parameter_ptr() << Core::IO::endl;
    }
  }

  // Conditions definition
  discret_ptr_->get_condition("LinePeriodicRve", line_periodic_rve_conditions_);
  discret_ptr_->get_condition("SurfacePeriodicRve", surface_periodic_rve_conditions_);
  discret_ptr_->get_condition("PointPeriodicRveReferenceNode", point_periodic_rve_ref_conditions_);
  discret_ptr_->get_condition(
      "PointLinearCoupledEquation", point_linear_coupled_equation_conditions_);

  // Input Checks: Dimensions
  if (line_periodic_rve_conditions_.size() == 0 && surface_periodic_rve_conditions_.size() != 0)
  {
    Core::IO::cout(Core::IO::verbose) << "Rve dimension: 3d" << Core::IO::endl;
    rve_dim_ = Inpar::RveMpc::RveDimension::rve3d;
  }
  else if (line_periodic_rve_conditions_.size() != 0 &&
           surface_periodic_rve_conditions_.size() == 0)
  {
    Core::IO::cout(Core::IO::verbose) << "Rve dimensions: 2d" << Core::IO::endl;
    rve_dim_ = Inpar::RveMpc::RveDimension::rve2d;
  }
  else
  {
    FOUR_C_THROW("Periodic rve edge condition cannot be combinded with peridodic rve surf cond. ");
  }

  // Input Checks
  Core::IO::cout(Core::IO::verbose)
      << "There are " << line_periodic_rve_conditions_.size()
      << " periodic rve edge conditions defined (2D)" << Core::IO::endl;
  if (line_periodic_rve_conditions_.size() != 0)
  {
    if (line_periodic_rve_conditions_.size() != 4 && line_periodic_rve_conditions_.size() != 2)
    {
      Core::IO::cout(Core::IO::verbose)
          << "Number of Conditions: " << line_periodic_rve_conditions_.size() << Core::IO::endl;
      FOUR_C_THROW("For a 2D RVE either all or two opposing edges must be used for PBCs");
    }
  }

  Core::IO::cout(Core::IO::verbose)
      << "There are " << surface_periodic_rve_conditions_.size()
      << " periodic rve surface conditions defined (3D)" << Core::IO::endl;

  Core::IO::cout(Core::IO::verbose)
      << "There are " << point_linear_coupled_equation_conditions_.size()
      << " linear coupled equations" << Core::IO::endl;

  if (point_periodic_rve_ref_conditions_.size() == 0 &&
      rve_ref_type_ == Inpar::RveMpc::RveReferenceDeformationDefinition::manual)
  {
    FOUR_C_THROW(
        "A DESIGN POINT PERIODIC RVE 2D BOUNDARY REFERENCE CONDITIONS is req. for manual ref. "
        "point "
        "definition");
  }

  if (point_periodic_rve_ref_conditions_.size() != 0 &&
      rve_ref_type_ == Inpar::RveMpc::RveReferenceDeformationDefinition::automatic)
    FOUR_C_THROW("Set the RVE_REFERENCE_POINTS to manual");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::find_opposite_edge_node(
    const int nodeID, Inpar::RveMpc::RveEdgeIdentifiers edge,
    std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap)
{
  std::string newPos;
  std::array<double, 2> R_ipim;

  Core::Nodes::Node* nodeA = discret_ptr_->g_node(nodeID);

  switch (edge)
  {
    case Inpar::RveMpc::RveEdgeIdentifiers::Gamma_xm:
    {
      R_ipim = r_xmxp_;
      newPos = "x+";
      Core::IO::cout(Core::IO::debug)
          << " Find partner node of Node " << nodeID << " on Edge:  " << newPos << Core::IO::endl;
      break;
    }

    case Inpar::RveMpc::RveEdgeIdentifiers::Gamma_ym:
    {
      R_ipim = r_ymyp_;
      newPos = "y+";
      Core::IO::cout(Core::IO::debug)
          << " Find partner node of Node " << nodeID << " on Edge:  " << newPos << Core::IO::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Specifiy the negative edge, 3D not implemented");
    }
  }

  // Calculate the Position of the opposing edge node
  std::vector<double> newPosition = {0.0, 0.0};
  for (int i = 0; i < 2; i++) newPosition[i] = nodeA->x()[i] + R_ipim[i];

  Core::IO::cout(Core::IO::debug) << "Calculated position of matching node: " << newPosition[0]
                                  << " / " << newPosition[1] << Core::IO::endl;

  // Loop all nodes of the relevant opposite edge line
  // ToDo: Switch to ArborX
  for (auto pairId : *rveBoundaryNodeIdMap[newPos])
  {
    Core::Nodes::Node* nodeB = discret_ptr_->g_node(pairId);

    if (std::abs(nodeB->x()[0] - newPosition[0]) < node_search_toler_)
    {
      if (std::abs(nodeB->x()[1] - newPosition[1]) < node_search_toler_)
      {
        Core::IO::cout(Core::IO::debug)
            << "Found Node Pair (IDs): A: " << nodeA->id() << " B: " << pairId << Core::IO::endl;
        return pairId;
      }
    }
  }
  FOUR_C_THROW("No matching node found! - Is the mesh perodic?");
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::build_periodic_mp_cs(
    std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap,
    std::map<std::string, int>& rveCornerNodeIdMap)
{
  std::vector<std::vector<Core::Nodes::Node*>> PBCs;
  std::vector<Core::Nodes::Node*> PBC;
  {
    Core::IO::cout(Core::IO::verbose)
        << "\nCreating Node Pairs for Periodic Boundary Conditions" << Core::IO::endl
        << "+--------------------------------------------------------------------+"
        << Core::IO::endl;
    Core::IO::cout(Core::IO::verbose) << "RVE Type: ";
  }
  switch (rve_ref_type_)
  {
    case Inpar::RveMpc::RveReferenceDeformationDefinition::automatic:
    {
      switch (rve_dim_)
      {
        case Inpar::RveMpc::RveDimension::rve3d:
        case Inpar::RveMpc::RveDimension::rve2d:
        {
          int numDim = 3;
          std::map<std::string, std::string> refEndNodeMap = {
              {"x", "N2"}, {"y", "N4"}, {"z", "N5"}};

          if (rve_dim_ == Inpar::RveMpc::rve3d)
          {
            Core::IO::cout(Core::IO::verbose) << "General 3D RVE" << Core::IO::endl;
          }
          else if (rve_dim_ == Inpar::RveMpc::rve2d)
          {
            Core::IO::cout(Core::IO::verbose) << "General 2D RVE" << Core::IO::endl;
            refEndNodeMap.erase("z");
            numDim = 2;
          }

          std::map<std::string, std::vector<double>> rveRefVecMap;
          std::vector<double> rveRefVector;

          for (const auto& surf : refEndNodeMap)
          {
            {
              rveRefVector.clear();
              for (int i = 0; i < numDim; ++i)
              {
                rveRefVector.push_back(
                    discret_ptr_->g_node(rveCornerNodeIdMap[surf.second])->x()[i] -
                    discret_ptr_->g_node(rveCornerNodeIdMap["N1"])->x()[i]);
              }
              rveRefVecMap[surf.first] = rveRefVector;

              Core::IO::cout(Core::IO::verbose)
                  << "RVE reference vector " << surf.first << " dimension: "
                  << "[ " << rveRefVector[0] << "; " << rveRefVector[1];
              if (rve_dim_ == Inpar::RveMpc::rve3d)
                Core::IO::cout(Core::IO::verbose) << "; " << rveRefVector[2];
              Core::IO::cout(Core::IO::verbose) << " ]" << Core::IO::endl;
            }
          }
          // Create PBC Node Pairs:
          // #ToDo: Use ArborX for the search
          std::string surfId;
          std::vector<double> matchPosition = {0.0, 0.0, 0.0};

          for (const auto& surf : refEndNodeMap)
          {
            for (auto nodeP : *rveBoundaryNodeIdMap[surf.first + "+"])
            {
              Core::IO::cout(Core::IO::debug)
                  << "Current SURF: " << surf.first << "Node+ ID: " << nodeP << Core::IO::endl;

              for (int i = 0; i < numDim; ++i)
              {
                matchPosition[i] =
                    discret_ptr_->g_node(nodeP)->x()[i] - rveRefVecMap[surf.first][i];
              }
              Core::IO::cout(Core::IO::debug)
                  << "Node+ location: " << discret_ptr_->g_node(nodeP)->x()[0] << ", "
                  << discret_ptr_->g_node(nodeP)->x()[1];
              if (rve_dim_ == Inpar::RveMpc::rve3d)
                Core::IO::cout(Core::IO::debug) << ", " << discret_ptr_->g_node(nodeP)->x()[2];
              Core::IO::cout(Core::IO::debug) << Core::IO::endl;

              Core::IO::cout(Core::IO::debug)
                  << "Position of matching Node: " << matchPosition[0] << ", " << matchPosition[1];
              if (rve_dim_ == Inpar::RveMpc::rve3d)
                Core::IO::cout(Core::IO::verbose) << ", " << matchPosition[2];
              Core::IO::cout(Core::IO::verbose) << Core::IO::endl;
              for (auto nodeM : *rveBoundaryNodeIdMap[surf.first + "-"])
              {
                if (nodeP == rveCornerNodeIdMap[surf.second]) break;

                bool isMatch = true;
                for (int i = 0; i < numDim; ++i)
                {
                  if (std::abs(discret_ptr_->g_node(nodeM)->x()[i] - matchPosition[i]) >
                      node_search_toler_)
                  {
                    isMatch = false;
                    break;
                  }
                }
                if (isMatch)
                {
                  PBC.push_back(discret_ptr_->g_node(nodeP));
                  PBC.push_back(discret_ptr_->g_node(nodeM));
                  PBC.push_back(discret_ptr_->g_node(rveCornerNodeIdMap[surf.second]));
                  PBC.push_back(discret_ptr_->g_node(rveCornerNodeIdMap["N1"]));
                  PBCs.push_back(PBC);
                  PBC.clear();

                  {
                    Core::IO::cout(Core::IO::debug)
                        << "nodes coupled on " << surf.first << " +/- boundary: " << nodeP << " - "
                        << nodeM << " = " << rveCornerNodeIdMap[surf.second] << " - "
                        << rveCornerNodeIdMap["N1"] << Core::IO::endl;
                  }
                }
              }
            }
          }
        }
        break;
      }
      break;
    }
    case Inpar::RveMpc::RveReferenceDeformationDefinition::manual:
    {
      Core::IO::cout(Core::IO::verbose) << "General 2D RVE" << Core::IO::endl;
      for (const auto& elem : rve_ref_node_map_)
      {
        Core::IO::cout(Core::IO::verbose) << elem.first << " first " << elem.second->id() << " sec "
                                          << "\n";
      }

      /* Loop over X- Edge */
      for (auto nodeXm : *rveBoundaryNodeIdMap["x-"])
      {
        if (nodeXm != rve_ref_node_map_["N1L"]->id())  // exlude N1 - N2 = N1 - N2
        {
          PBC.push_back(discret_ptr_->g_node(nodeXm));
          PBC.push_back(discret_ptr_->g_node(find_opposite_edge_node(
              nodeXm, Inpar::RveMpc::RveEdgeIdentifiers::Gamma_xm, rveBoundaryNodeIdMap)));
          PBC.push_back(rve_ref_node_map_["N1L"]);
          PBC.push_back(rve_ref_node_map_["N2"]);
          PBCs.push_back(PBC);

          Core::IO::cout(Core::IO::debug)
              << "PBC MPC Set created X- ---> X+ Edge:  " << Core::IO::endl;
          for (auto* nnn : PBC)
          {
            Core::IO::cout(Core::IO::debug) << " ___ " << nnn->id();
          }
          Core::IO::cout(Core::IO::debug) << Core::IO::endl;
          PBC.clear();
        }
      }

      /* Loop over Y- Edge*/
      for (auto nodeYm : *rveBoundaryNodeIdMap["y-"])
      {
        if (nodeYm != rve_ref_node_map_["N1B"]->id() && nodeYm != rve_ref_node_map_["N2"]->id())
        {
          PBC.push_back(discret_ptr_->g_node(find_opposite_edge_node(
              nodeYm, Inpar::RveMpc::RveEdgeIdentifiers::Gamma_ym, rveBoundaryNodeIdMap)));
          PBC.push_back(discret_ptr_->g_node(nodeYm));
          PBC.push_back(rve_ref_node_map_["N4"]);
          PBC.push_back(rve_ref_node_map_["N1B"]);
          PBCs.push_back(PBC);
          PBC.clear();
        }
      }
      break;
    }
    default:
      FOUR_C_THROW("No ref def type defined");
  }

  // Ensure no constraint is enforced twice:
  int indx = 0;
  std::vector<int> ids;
  std::map<int, std::vector<int>> idListSet;

  for (const auto& pbc : PBCs)
  {
    for (auto* node : pbc)
    {
      ids.push_back(node->id());
    }
    std::sort(ids.begin(), ids.end());
    idListSet[indx++] = ids;
    ids.clear();
  }


  // Print
  int nr = 0;
  Core::IO::cout(Core::IO::debug) << "Sorted Pair Ids:" << Core::IO::endl;
  for (const auto& couple : idListSet)
  {
    Core::IO::cout(Core::IO::debug) << "Pair " << nr++ << ":";
    for (auto id : couple.second)
    {
      Core::IO::cout(Core::IO::debug) << id << ", ";
    }
    Core::IO::cout(Core::IO::debug) << Core::IO::endl;
  }

  std::vector<int> idsToRemove;
  for (const auto& entryA : idListSet)
  {
    if (std::find(idsToRemove.begin(), idsToRemove.end(), entryA.first) == idsToRemove.end())
    {
      for (const auto& entryB : idListSet)
      {
        if (entryA.second == entryB.second && entryA.first != entryB.first)
        {
          idsToRemove.push_back(entryB.first);
          Core::IO::cout(Core::IO::debug) << "remove pair: " << entryB.first << Core::IO::endl;
        }
      }
    }
  }

  // Remove duplicate constraints
  std::sort(idsToRemove.rbegin(), idsToRemove.rend());
  for (int id : idsToRemove)
  {
    PBCs.erase(PBCs.begin() + id);
  }
  Core::IO::cout(Core::IO::verbose)
      << "All Node Pairs found. Following Nodes are coupled:" << Core::IO::endl;

  // Create the vector of MPC "Elements
  int mpcId = 0;
  std::vector<int> pbcDofs;
  std::vector<double> pbcCoefs = {1., -1., -1., 1.};

  int nDofCoupled = 3;
  if (rve_dim_ == Inpar::RveMpc::rve2d)
  {
    nDofCoupled = 2;
  }


  for (const auto& pbc : PBCs)
  {
    Core::IO::cout(Core::IO::debug) << "Node Pair: ";
    for (int dim = 0; dim < nDofCoupled; ++dim)
    {
      pbcDofs.clear();
      for (auto* node : pbc)
      {
        Core::IO::cout(Core::IO::debug) << node->id() << ", ";

        // Create coupled equation dof list:
        pbcDofs.emplace_back(discret_ptr_->dof(node)[dim]);
      }
      Core::IO::cout(Core::IO::debug) << "\ndofs used in linear coupled equation: \n ";
      for (auto dof : pbcDofs)
      {
        Core::IO::cout(Core::IO::debug) << dof << ", ";
      }
      Core::IO::cout(Core::IO::debug) << Core::IO::endl;
      listMPCs_.emplace_back(std::make_shared<LinearCoupledEquation>(mpcId++, pbcDofs, pbcCoefs));
    }
    Core::IO::cout(Core::IO::debug) << "\n";
  }
  Core::IO::cout(Core::IO::verbose)
      << Core::IO::endl
      << "Total number of node pairs created for periodic boundary conditions: " << listMPCs_.size()
      << Core::IO::endl;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::build_linear_mp_cs()
{
  Core::IO::cout(Core::IO::verbose)
      << Core::IO::endl
      << "-------------------------------------------" << Core::IO::endl;
  Core::IO::cout(Core::IO::verbose)
      << "Reading linear coupled eq. from .dat file" << Core::IO::endl;
  Core::IO::cout(Core::IO::verbose)
      << "linear MPC codition count: " << point_linear_coupled_equation_conditions_.size()
      << Core::IO::endl;


  int nEq = 0;
  for (const auto& ceTerm : point_linear_coupled_equation_conditions_)
  {
    nEq = std::max(nEq, (ceTerm->parameters().get<int>("EQUATION")) - 1);
  }
  Core::IO::cout(Core::IO::verbose)
      << "There are " << nEq + 1 << " linear MPC Equations defined" << Core::IO::endl;

  int dofPos;
  int cond_num = 0;
  std::vector<std::vector<int>> constraintRowIds(nEq + 1);
  std::vector<std::vector<int>> constraintColIds(nEq + 1);
  std::vector<std::vector<double>> constraintCoeffs(nEq + 1);


  for (const auto& ceTerm : point_linear_coupled_equation_conditions_)
  {
    auto eq_id = (ceTerm->parameters().get<int>("EQUATION")) - 1;
    const auto* node_id = ceTerm->get_nodes();
    const auto& dofStr = ceTerm->parameters().get<std::string>("ADD");
    auto coef = ceTerm->parameters().get<double>("COEFFICIENT");
    auto* node = discret_ptr_->g_node(node_id->data()[0]);


    if (dofStr == "dispx")
    {
      dofPos = 0;
    }
    else if (dofStr == "dispy")
    {
      dofPos = 1;
    }
    else
    {
      FOUR_C_THROW(
          "Linear coupled equations (MPCs) are only implemented for 2D (dispx or dispy DOFs)");
    }
    auto dofID = discret_ptr_->dof(node)[dofPos];


    Core::IO::cout(Core::IO::debug) << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "Condition Number " << cond_num++ << ": " << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "Eq.Id: " << eq_id << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "Node Id: " << node_id->data()[0] << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "Disp String: " << dofStr.c_str() << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "DOF ID: " << dofID << Core::IO::endl;
    Core::IO::cout(Core::IO::debug) << "COEF: " << coef << Core::IO::endl << Core::IO::endl;

    // Save the linear MPCs
    constraintRowIds[eq_id].push_back(eq_id);
    constraintColIds[eq_id].push_back(dofID);
    constraintCoeffs[eq_id].push_back(coef);
    Core::IO::cout(Core::IO::debug) << "Added Term Equation with ID: " << eq_id << Core::IO::endl;
    Core::IO::cout(Core::IO::debug)
        << "Current SIze constraintColIDs" << constraintColIds.size() << Core::IO::endl;
  }
  // Get number of MPC already in the MPCs List
  int nMPC = 0;
  for (const auto& mpc : listMPCs_)
  {
    nMPC += mpc->get_number_of_mp_cs();
  }
  unsigned int i = 0;
  for (; i < constraintRowIds.size(); ++i)
  {
    listMPCs_.emplace_back(
        std::make_shared<LinearCoupledEquation>(nMPC++, constraintColIds[i], constraintCoeffs[i]));

    Core::IO::cout(Core::IO::verbose) << "Linear MPC #" << i << "  Created: 0 = ";
    for (unsigned int o = 0; o < constraintColIds[i].size(); ++o)
    {
      Core::IO::cout(Core::IO::verbose)
          << " +" << constraintCoeffs[i][o] << "*d" << constraintColIds[i][o];
    }
    Core::IO::cout(Core::IO::verbose) << Core::IO::endl;
  }
  Core::IO::cout(Core::IO::verbose) << "Number of Linear MPCs Created: " << i + 1 << Core::IO::endl;
  Core::IO::cout(Core::IO::verbose)
      << "Number of Elements in the listMPC: " << nMPC << Core::IO::endl;

  return i + 1;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::find_periodic_rve_corner_nodes(
    const std::vector<int>* edge1, const std::vector<int>* edge2)
{
  for (int nodeId : *edge2)
  {
    if (std::find(edge1->begin(), edge1->end(), nodeId) != edge1->end()) return nodeId;
  }
  return -1;
}

int CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::find_periodic_rve_corner_nodes(
    const std::vector<int>* surf1, const std::vector<int>* surf2, const std::vector<int>* surf3)
{
  std::vector<int> commonNodeIds12, commonNodeIds123;
  std::set_intersection(surf1->begin(), surf1->end(), surf2->begin(), surf2->end(),
      std::back_inserter(commonNodeIds12));

  std::set_intersection(surf3->begin(), surf3->end(), commonNodeIds12.begin(),
      commonNodeIds12.end(), std::back_inserter(commonNodeIds123));

  if (commonNodeIds123.size() < 1)
  {
    FOUR_C_THROW("No common node found");
  }

  else if (commonNodeIds123.size() > 1)
  {
    FOUR_C_THROW("More than one common node found");
  }
  return commonNodeIds123[0];
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::
    build_periodic_rve_corner_node_map(
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap,
        std::map<std::string, int>& rveCornerNodeIdMap)
{
  switch (rve_dim_)
  {
    case Inpar::RveMpc::RveDimension::rve2d:
    {
      //* Get the Corner Node Ids */
      /*              N4 -- N
       *              |     |
       *      N   --- |     |---   N
       *      |                    |
       *    x-|                    |x+
       *      N1L ---|       | ----N2
       *             |       |
       *             N1B --  N
       *
       *      N1 ------ y- ------ N2
       */

      //* Get the Corner Node Ids */

      /*
       *      N4 ------y+------ N3
       *      |                    |
       *      x-                   |x+
       *      |                   |
       *      N1 ------ y- ------ N2
       */

      rveCornerNodeIdMap["N1"] =
          find_periodic_rve_corner_nodes(rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["y-"]);
      rveCornerNodeIdMap["N2"] =
          find_periodic_rve_corner_nodes(rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["y-"]);
      rveCornerNodeIdMap["N3"] =
          find_periodic_rve_corner_nodes(rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["y+"]);
      rveCornerNodeIdMap["N4"] =
          find_periodic_rve_corner_nodes(rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["y+"]);


      Core::IO::cout(Core::IO::verbose)
          << "\nAutomatically determined following reference Nodes: " << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose)
          << "---------------------------------------------------" << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << "N1: " << rveCornerNodeIdMap["N1"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N2: " << rveCornerNodeIdMap["N2"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N3: " << rveCornerNodeIdMap["N3"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N4: " << rveCornerNodeIdMap["N4"] << Core::IO::endl;
    }
    break;
    case Inpar::RveMpc::RveDimension::rve3d:
    {
      //  z ^        N8 +  +   +   N7
      //    |      + .            + +
      //    |    +   .          +   +
      //    |  +     .        +     +
      //    N5 +  +  +  +  N3       +
      //    +       .      +        +
      //    +       N4 .  .+   . .. N3
      //    +      .       +     +
      //    +    .         +  +
      //    +  .           +
      //   N1 +  +  +  +  N2 ------> x
      //

      std::vector<std::string> boundaryNames = {"x+", "x-", "y+", "y-", "z+", "z-"};
      rveCornerNodeIdMap["N1"] = find_periodic_rve_corner_nodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y-"]);

      rveCornerNodeIdMap["N2"] = find_periodic_rve_corner_nodes(
          rveBoundaryNodeIdMap["x+"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y-"]);

      rveCornerNodeIdMap["N4"] = find_periodic_rve_corner_nodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z-"], rveBoundaryNodeIdMap["y+"]);

      rveCornerNodeIdMap["N5"] = find_periodic_rve_corner_nodes(
          rveBoundaryNodeIdMap["x-"], rveBoundaryNodeIdMap["z+"], rveBoundaryNodeIdMap["y-"]);

      Core::IO::cout(Core::IO::verbose)
          << "\nAutomatically determined following reference Nodes: " << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose)
          << "---------------------------------------------------" << Core::IO::endl;
      Core::IO::cout(Core::IO::verbose) << "N1: " << rveCornerNodeIdMap["N1"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N2: " << rveCornerNodeIdMap["N2"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N4: " << rveCornerNodeIdMap["N4"] << ", ";
      Core::IO::cout(Core::IO::verbose) << "N5: " << rveCornerNodeIdMap["N5"] << Core::IO::endl;
      break;
    }
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::
    build_periodic_rve_boundary_node_map(
        std::map<std::string, const std::vector<int>*>& rveBoundaryNodeIdMap)
{
  switch (rve_dim_)
  {
    case Inpar::RveMpc::RveDimension::rve2d:
    {
      discret_ptr_->get_condition("LinePeriodicRve", line_periodic_rve_conditions_);



      Core::IO::cout(Core::IO::verbose)
          << Core::IO::endl
          << "Reading Line Conditions:  " << Core::IO::endl
          << "+--------------------------------------------------------------------+"
          << Core::IO::endl;
      for (const auto& conditionLine : line_periodic_rve_conditions_)
      {
        const auto& boundary = conditionLine->parameters().get<std::string>("EDGE");

        // Print the Edge Condition
        Core::IO::cout(Core::IO::verbose) << "EDGE: " << boundary.c_str() << " Node IDs: ";
        for (auto nodeId : *conditionLine->get_nodes())
        {
          Core::IO::cout(Core::IO::verbose) << nodeId << " ";
        }
        Core::IO::cout(Core::IO::verbose) << Core::IO::endl;

        // Create EdgeNodeMap
        rveBoundaryNodeIdMap[boundary.c_str()] = conditionLine->get_nodes();
      }
    }
    break;

    case Inpar::RveMpc::RveDimension::rve3d:
    {
      discret_ptr_->get_condition("SurfacePeriodicRve", surface_periodic_rve_conditions_);

      Core::IO::cout(Core::IO::verbose)
          << Core::IO::endl
          << "Reading Surface Conditions:  " << Core::IO::endl
          << "+--------------------------------------------------------------------+"
          << Core::IO::endl;
      for (const auto& conditionLine : surface_periodic_rve_conditions_)
      {
        const auto& boundary = conditionLine->parameters().get<std::string>("SURF");

        Core::IO::cout(Core::IO::verbose) << "SURF: " << boundary.c_str() << " Node IDs: ";
        for (auto nodeId : *conditionLine->get_nodes())
        {
          Core::IO::cout(Core::IO::verbose) << nodeId << " ";
        }
        Core::IO::cout(Core::IO::verbose) << Core::IO::endl;
        rveBoundaryNodeIdMap[boundary.c_str()] = conditionLine->get_nodes();
      }
      break;
    }
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::reset() {}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

std::map<Solid::EnergyType, double>
CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager::get_energy() const
{
  FOUR_C_THROW("This function is not implemented for the RveMultiPointConstraintManager.");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE
