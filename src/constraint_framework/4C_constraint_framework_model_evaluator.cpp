// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_constraint_framework_model_evaluator.hpp"

#include "4C_constraint_framework_submodelevaluator_embeddedmesh.hpp"
#include "4C_constraint_framework_submodelevaluator_mpc.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::setup()
{
  check_init();

  constraint_stiff_ptr_ =
      std::make_shared<Core::LinAlg::SparseMatrix>(*global_state().get_discret()->dof_row_map(), 81,
          true, true, Core::LinAlg::SparseMatrix::FE_MATRIX);

  constraint_force_ptr_ =
      std::make_shared<Core::LinAlg::Vector<double>>(*global_state().dof_row_map_view(), true);

  set_sub_model_types();
  create_sub_model_evaluators();

  visualization_params_ = Core::IO::visualization_parameters_factory(
      Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
      *Global::Problem::instance()->output_control_file(), global_state().get_time_n());

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::set_sub_model_types()
{
  check_init();

  submodeltypes_ = std::set<enum Inpar::CONSTRAINTS::SubModelType>();

  // ---------------------------------------------------------------------------
  // check for multi point constraints
  // ---------------------------------------------------------------------------
  std::vector<std::shared_ptr<Core::Conditions::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation, embeddedMeshConditions;

  discret_ptr()->get_condition("LinePeriodicRve", linePeriodicRve);
  discret_ptr()->get_condition("SurfacePeriodicRve", surfPeriodicRve);
  discret_ptr()->get_condition("PointLinearCoupledEquation", pointLinearCoupledEquation);
  discret_ptr()->get_condition("EmbeddedMeshSolidSurfCoupling", embeddedMeshConditions);

  if (linePeriodicRve.size() > 0 || surfPeriodicRve.size() > 0 ||
      pointLinearCoupledEquation.size() > 0)
  {
    submodeltypes_.insert(Inpar::CONSTRAINTS::SubModelType::submodel_pbc_rve);
  }
  if (embeddedMeshConditions.size() > 0)
  {
    submodeltypes_.insert(Inpar::CONSTRAINTS::SubModelType::submodel_embeddedmesh);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::create_sub_model_evaluators()
{
  // Create vector with the Sub-model-evaluators
  sub_model_vec_ptr_ = Solid::ModelEvaluator::Constraints::SubmodelevaluatorVector(0);

  for (const auto& submodeltype : submodeltypes_)
  {
    switch (submodeltype)
    {
      case Inpar::CONSTRAINTS::SubModelType::submodel_pbc_rve:
      {
        sub_model_vec_ptr_.emplace_back(
            std::make_shared<CONSTRAINTS::SUBMODELEVALUATOR::RveMultiPointConstraintManager>(
                discret_ptr(), constraint_stiff_ptr_.get()));

        break;
      }
      case Inpar::CONSTRAINTS::SubModelType::submodel_embeddedmesh:
      {
        sub_model_vec_ptr_.emplace_back(
            std::make_shared<CONSTRAINTS::SUBMODELEVALUATOR::EmbeddedMeshConstraintManager>(
                discret_ptr(), *global_state().get_dis_np().get()));

        break;
      }
      default:
      {
        FOUR_C_THROW(
            "Something went wrong: Apparently a Constraint ME was created that is not "
            "required. Check the Adapter");
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::have_sub_model_type(
    Inpar::CONSTRAINTS::SubModelType const& submodeltype) const
{
  check_init();
  return (submodeltypes_.find(submodeltype) != submodeltypes_.end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::reset(const Core::LinAlg::Vector<double>& x)
{
  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->reset();
  }
  constraint_stiff_ptr_->zero();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::evaluate_force()
{
  pre_evaluate();
  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->evaluate_force_stiff(
        *global_state().get_dis_np().get(), global_state_ptr(), nullptr, constraint_force_ptr_);
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::evaluate_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->un_complete();
  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->evaluate_force_stiff(
        *global_state().get_dis_np().get(), global_state_ptr(), constraint_stiff_ptr_, nullptr);
  }
  if (not constraint_stiff_ptr_->filled()) constraint_stiff_ptr_->complete();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::evaluate_force_stiff()
{
  pre_evaluate();

  constraint_stiff_ptr_->un_complete();
  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->evaluate_force_stiff(*global_state().get_dis_np().get(), global_state_ptr(),
        constraint_stiff_ptr_, constraint_force_ptr_);
  }
  if (not constraint_stiff_ptr_->filled()) constraint_stiff_ptr_->complete();

  return true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::pre_evaluate()
{
  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->evaluate_coupling_terms(*global_state_ptr());
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::assemble_force(
    Core::LinAlg::Vector<double>& f, const double& timefac_np) const
{
  Core::LinAlg::assemble_my_vector(1.0, f, timefac_np, *constraint_force_ptr_);
  constraint_force_ptr_->PutScalar(0.0);
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Solid::ModelEvaluator::Constraints::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  std::shared_ptr<Core::LinAlg::SparseMatrix> jac_dd_ptr = global_state().extract_displ_block(jac);

  jac_dd_ptr->add(*constraint_stiff_ptr_, false, timefac_np, 1.0);

  constraint_stiff_ptr_->zero();
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // There is nothing to write for now
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  // There is nothing to read for now
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::predict(const Inpar::Solid::PredEnum& pred_type) {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::update_step_state(const double& timefac_n)
{
  if (constraint_force_ptr_)
  {
    std::shared_ptr<Core::LinAlg::Vector<double>>& fstruct_ptr =
        global_state().get_fstructure_old();
    fstruct_ptr->Update(timefac_n, *constraint_force_ptr_, 1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::update_step_element() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_stress_strain() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_energy()
{
  check_init_setup();

  std::map<Solid::EnergyType, double> energy_this_submodel;

  for (auto& sme_iter : sub_model_vec_ptr_)
  {
    energy_this_submodel = sme_iter->get_energy();

    for (auto const& energy_type : energy_this_submodel)
      eval_data().add_contribution_to_energy_type(energy_type.second, energy_type.first);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::determine_optional_quantity() {}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::reset_step_state()
{
  FOUR_C_THROW("This function is not implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::runtime_pre_output_step_state() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::runtime_output_step_state() const
{
  // Write output vtk of lagrange multipliers
  std::pair<double, int> output_time_and_step;
  if (visualization_params_.every_iteration_ == true)
  {
    output_time_and_step = Core::IO::get_time_and_time_step_index_for_output(visualization_params_,
        global_state().get_time_n(), global_state().get_step_n(), eval_data().get_nln_iter());
  }
  else
  {
    output_time_and_step = Core::IO::get_time_and_time_step_index_for_output(
        visualization_params_, global_state().get_time_n(), global_state().get_step_n());
  }

  for (auto& some_iter : sub_model_vec_ptr_)
  {
    some_iter->runtime_output_step_state(output_time_and_step);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<const Epetra_Map> Solid::ModelEvaluator::Constraints::get_block_dof_row_map_ptr()
    const
{
  return global_state().dof_row_map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Solid::ModelEvaluator::Constraints::get_current_solution_ptr() const
{
  // there are no model specific solution entries
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Solid::ModelEvaluator::Constraints::get_last_time_step_solution_ptr() const
{
  // there are no model specific solution entries
  return nullptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::post_output() {}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::
    evaluate_jacobian_contributions_from_element_level_for_ptc()
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::assemble_jacobian_contributions_from_element_level_for_ptc(
    std::shared_ptr<Core::LinAlg::SparseMatrix>& modjac, const double& timefac_n)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::create_backup_state(
    const Core::LinAlg::Vector<double>& dir)
{
  FOUR_C_THROW("This function is not yet implemented");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::ModelEvaluator::Constraints::recover_from_backup_state()
{
  FOUR_C_THROW("This function is not yet implemented");
}

FOUR_C_NAMESPACE_CLOSE
