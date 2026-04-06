// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/invm_manager.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/mf_operator.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>

#include <prismspf/config.h>

#include <memory>
//
#include <deal.II/lac/precondition_block.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolveContext;

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree, typename number>
class LinearSolver : public SolverBase<dim, degree, number>
{
protected:
  using SolverBase<dim, degree, number>::solutions;
  using SolverBase<dim, degree, number>::solve_context;
  using SolverBase<dim, degree, number>::solve_block;

public:
  /**
   * @brief Constructor.
   */
  LinearSolver(SolveBlock                               _solve_block,
               const SolveContext<dim, degree, number> &_solve_context)
    : SolverBase<dim, degree, number>(_solve_block, _solve_context)
    , lin_params(
        solve_context->get_user_inputs().linear_solve_parameters.linear_solvers.at(
          solve_block.id))
  {}

  /**
   * @brief Initialize the solver.
   */
  void
  init(const std::list<DependencyMap> &all_dependeny_sets) override
  {
    SolverBase<dim, degree, number>::init(all_dependeny_sets);
    unsigned int num_levels = solve_context->get_dof_manager().get_dof_handlers().size();
    rhs_vector.resize(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_vector[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
    // Initialize rhs_operators
    rhs_operators.reserve(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_operators.emplace_back(solve_context->get_pde_operator(),
                                   &PDEOperatorBase<dim, degree, number>::compute_rhs,
                                   solve_context->get_field_attributes(),
                                   solve_context->get_solution_indexer(),
                                   relative_level,
                                   solve_block.dependencies_rhs,
                                   solve_context->get_simulation_timer());
        rhs_operators[relative_level].initialize(solutions);
        rhs_operators[relative_level].set_scaling_diagonal(
          lin_params.tolerance_type != AbsoluteResidual,
          solve_context->get_invm_manager().get_invm_sqrt(
            solve_context->get_field_attributes(),
            solve_block.field_indices,
            relative_level));
      }
    // Initialize lhs_operators
    lhs_operators.reserve(num_levels);
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        lhs_operators.emplace_back(solve_context->get_pde_operator(),
                                   &PDEOperatorBase<dim, degree, number>::compute_lhs,
                                   solve_context->get_field_attributes(),
                                   solve_context->get_solution_indexer(),
                                   relative_level,
                                   solve_block.dependencies_lhs,
                                   solve_context->get_simulation_timer());
        lhs_operators[relative_level].initialize(solutions);
        lhs_operators[relative_level].set_scaling_diagonal(
          lin_params.tolerance_type != AbsoluteResidual,
          solve_context->get_invm_manager().get_invm_sqrt(
            solve_context->get_field_attributes(),
            solve_block.field_indices,
            relative_level));
      }
    linear_solver_control.set_max_steps(lin_params.max_iterations);
    linear_solver_control.set_tolerance(lin_params.tolerance * normalization_value());
    inhomogenous_values.reinit(solutions.get_solution_full_vector(0));
    solutions.apply_constraints(inhomogenous_values, 0);
    inhomogenous_rhs.reinit(solutions.get_solution_full_vector(0));
  }

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    SolverBase<dim, degree, number>::reinit();
    const unsigned int num_levels = rhs_vector.size();
    for (unsigned int relative_level = 0; relative_level < num_levels; ++relative_level)
      {
        rhs_vector[relative_level].reinit(
          solutions.get_solution_full_vector(relative_level));
      }
    inhomogenous_values.reinit(solutions.get_solution_full_vector(0));
    solutions.apply_constraints(inhomogenous_values, 0);
    inhomogenous_rhs.reinit(solutions.get_solution_full_vector(0));
  }

  /**
   * @brief Solve for a single update step.
   */
  void
  solve_level(unsigned int relative_level) override
  {
    // Zero out the ghosts
    Timer::start_section("Zero ghosts");
    solutions.zero_out_ghosts(relative_level);
    Timer::end_section("Zero ghosts");

    // Set up rhs vector
    rhs_operators[relative_level].compute_operator(rhs_vector[relative_level]);
    if (relative_level == 0)
      {
        lhs_operators[0].read_plain = true;
        lhs_operators[0].compute_operator(inhomogenous_rhs, inhomogenous_values);
        lhs_operators[0].read_plain = false;
        rhs_vector[0] -= inhomogenous_rhs;
      }
    // Linear solve
    do_linear_solve(rhs_vector[relative_level],
                    lhs_operators[relative_level],
                    solutions.get_solution_full_vector(relative_level));

    if (relative_level == 0)
      {
        solutions.get_solution_full_vector(0) += inhomogenous_values;
      }
    // Apply constraints
    solutions.apply_constraints(relative_level);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    solutions.update_ghosts(relative_level);
    Timer::end_section("Update ghosts");
  }

  int
  do_linear_solve(BlockVector<number>             &b_vector,
                  MFOperator<dim, degree, number> &lhs_operator,
                  BlockVector<number>             &x_vector)
  {
    // Linear solve
    try
      {
        dealii::SolverCG<BlockVector<number>> cg_solver(linear_solver_control);
        if (lin_params.preconditioner == PreconditionerType::GMG)
          {
            cg_solver.solve(lhs_operator, x_vector, b_vector, *multigrid_preconditioner);
          }
        else
          {
            cg_solver.solve(lhs_operator,
                            x_vector,
                            b_vector,
                            dealii::PreconditionIdentity());
          }
        if (solve_context->get_user_inputs().output_parameters.should_output(
              solve_context->get_simulation_timer().get_increment()))
          {
            ConditionalOStreams::pout_summary()
              << " Linear solve final residual : "
              << linear_solver_control.last_value() / normalization_value()
              << " Linear steps: " << linear_solver_control.last_step() << "\n"
              << std::flush;
          }
      }
    catch (...) // TODO: more specific catch
      {
        ConditionalOStreams::pout_base()
          << "[Increment " << solve_context->get_simulation_timer().get_increment()
          << "] "
          << "Warning: linear solver did not converge as per set tolerances before "
          << lin_params.max_iterations << " iterations.\n";
      }
    return linear_solver_control.last_step();
  }

protected:
  /**
   * @brief Matrix free operators for each level
   */
  std::vector<MFOperator<dim, degree, number>> rhs_operators;

  std::vector<MFOperator<dim, degree, number>> lhs_operators;
  std::vector<BlockVector<number>>             rhs_vector;

  double
  normalization_value()
  {
    SolverToleranceType type  = lin_params.tolerance_type;
    double              value = 1.0;
    if (type == RMSEPerField || type == RMSETotal)
      {
        value *= std::sqrt(solve_context->get_triangulation_manager().get_volume());
      }
    if (type == RMSEPerField || type == IntegratedPerField)
      {
        value *= std::sqrt(double(solve_block.field_indices.size()));
      }
    return value;
  }

private:
  using MGTransferType =
    dealii::MGTransferBlockGlobalCoarsening<dim, BlockVector<number>>;
  /**
   * @brief Solver control. Contains max iterations and tolerance.
   */
  dealii::SolverControl linear_solver_control;
  /**
   * @brief Linear solve settings
   */
  LinearSolverParameters lin_params;

  /**
   * @brief Vector containing only the inhomogeneous constraints (namely, non-zero
   * Dirichlet values)
   */
  BlockVector<number> inhomogenous_values;

  /**
   * @brief Result of the linear operator applied to the inhomogeneous values.
   */
  BlockVector<number> inhomogenous_rhs;

  /**
   * @brief Multigrid preconditioner
   */
  std::shared_ptr<dealii::PreconditionMG<dim, BlockVector<number>, MGTransferType>>
    multigrid_preconditioner;

  void
  init_multigrid()
  {
    const unsigned int min_level = lin_params.min_mg_level;
    const unsigned int max_level =
      solve_context->get_user_inputs().spatial_discretization.max_refinement;

    // 1. Level operators
    dealii::MGLevelObject<MFOperator<dim, degree, number>> mg_lhs_operators(
      min_level,
      max_level,
      solve_context->get_pde_operator(),
      &PDEOperatorBase<dim, degree, number>::compute_lhs,
      solve_context->get_field_attributes(),
      solve_context->get_solution_indexer(),
      0,
      solve_block.dependencies_lhs,
      solve_context->get_simulation_timer());
    for (unsigned level = min_level; level < max_level; ++level)
      {
        const unsigned int relative_level = max_level - level;
        mg_lhs_operators[level]           = lhs_operators[relative_level];
        /* mg_lhs_operators[level]           = MFOperator<dim, degree, number>(
          solve_context->get_pde_operator(),
          &PDEOperatorBase<dim, degree, number>::compute_lhs,
          solve_context->get_field_attributes(),
          solve_context->get_solution_indexer(),
          relative_level,
          solve_block.dependencies_lhs,
          solve_context->get_simulation_timer());
        mg_lhs_operators[level].initialize(solutions);
        mg_lhs_operators[level].set_scaling_diagonal(
          true,
          solve_context->get_invm_manager().get_invm_sqrt(
            solve_context->get_field_attributes(),
            solve_block.field_indices,
            relative_level)); */
      }

    // 3. MG transfer
    // dealii::MGTransferGlobalCoarsening<dim, BlockVector<number>> ?
    // MGTransferBlockGlobalCoarsening ?
    // MGTransferBlockMatrixFree ?

    dealii::MGTransferMF<dim, number> mg_trans_mf;
    MGTransferType                    mg_transfer(mg_trans_mf); // Constraints?
    // NOTE: dof_handler.distribute_mg_dofs() must have been called
    mg_transfer.build(
      solve_context->get_dof_manager().get_field_dof_handlers(solve_block.field_indices,
                                                              0));

    // 4. MG Smoother (takes in operators) This is similar to a solver, but is
    // conceptually different.
    // Preconditioner for smoother.
    // TODO: use PreconditionBlockJacobi or other block preconditioner instead of
    // PreconditionChebyshev
    using SmootherPrecond =
      dealii::PreconditionChebyshev<MFOperator<dim, degree, number>, BlockVector<number>>;
    dealii::MGLevelObject<typename SmootherPrecond::AdditionalData> smoother_data(
      min_level,
      max_level);
    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        smoother_data[level].smoothing_range     = lin_params.smoothing_range;
        smoother_data[level].degree              = lin_params.smoother_degree;
        smoother_data[level].eig_cg_n_iterations = lin_params.eig_cg_n_iterations;
        // smoother_data[level].preconditioner; // todo
        smoother_data[level].constraints.close();

        unsigned int        relative_level = level - min_level;
        BlockVector<number> identity;
        identity.reinit(solutions.get_solution_full_vector(relative_level));
        for (int index : identity.locally_owned_elements())
          {
            identity[index] = 1.0;
          }
        smoother_data[level].preconditioner =
          std::make_shared<dealii::DiagonalMatrix<BlockVector<number>>>(identity);
      }
    // Wrapper around a generic preconditioner to be used as a smoother
    using Smoother = dealii::MGSmootherPrecondition<MFOperator<dim, degree, number>,
                                                    SmootherPrecond,
                                                    BlockVector<number>>;
    Smoother mg_smoother;
    mg_smoother.initialize(mg_lhs_operators, smoother_data);

    // 5. Coarse grid solver
    dealii::MGCoarseGridApplySmoother<BlockVector<number>> mg_coarse_solver;
    mg_coarse_solver.initialize(mg_smoother);

    // 6. Multigrid object
    dealii::mg::Matrix<BlockVector<number>> mg_matrix(mg_lhs_operators);
    dealii::Multigrid<BlockVector<number>>  multigrid(
      mg_matrix,
      mg_coarse_solver,
      mg_transfer,
      mg_smoother,
      mg_smoother,
      min_level,
      max_level,
      dealii::Multigrid<BlockVector<number>>::Cycle::v_cycle);

    // 7. Turn MG into a preconditioner object
    multigrid_preconditioner =
      std::make_shared<dealii::PreconditionMG<dim, BlockVector<number>, MGTransferType>>(
        solve_context->get_dof_manager().get_field_dof_handlers(solve_block.field_indices,
                                                                0),
        multigrid,
        mg_transfer);
  }
};

PRISMS_PF_END_NAMESPACE
