// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>

#include <prismspf/core/simulation_time.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, unsigned int degree, typename number>
class VariableContainer;

/**
 * @brief This class contains the user implementation of each PDE operator.
 */
template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator
{
public:
  using SizeType = dealii::VectorizedArray<number>;

  /**
   * @brief Constructor.
   */
  explicit PDEOperator(const UserInputParameters<dim> &user_inputs)
    : delta_t(user_inputs.get_temporal_discretization().get_timestep())
  {}

  /**
   * @brief Destructor.
   */
  virtual ~PDEOperator() = default;

  /**
   * @brief User-implemented class for the setting initial conditions.
   */
  virtual void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        number                   &scalar_value,
                        number                   &vector_component_value) const = 0;

  /**
   * @brief User-implemented class for the setting nonuniform boundary conditions.
   */
  virtual void
  set_nonuniform_dirichlet(const unsigned int       &index,
                           const unsigned int       &boundary_id,
                           const unsigned int       &component,
                           const dealii::Point<dim> &point,
                           number                   &scalar_value,
                           number                   &vector_component_value) const = 0;

  /**
   * @brief User-implemented class for the RHS of explicit equations.
   */
  virtual void
  compute_explicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                       const SimulationTime                   &simulation_time,
                       const dealii::Point<dim, SizeType>     &q_point_loc,
                       const SizeType                         &element_volume,
                       Types::Index                            solve_block) const = 0;

  /**
   * @brief User-implemented class for the RHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                          const SimulationTime                   &simulation_time,
                          const dealii::Point<dim, SizeType>     &q_point_loc,
                          const SizeType                         &element_volume,
                          Types::Index                            solve_block,
                          Types::Index index = Numbers::invalid_index) const = 0;

  /**
   * @brief User-implemented class for the LHS of nonexplicit equations.
   */
  virtual void
  compute_nonexplicit_lhs(VariableContainer<dim, degree, number> &variable_list,
                          const SimulationTime                   &simulation_time,
                          const dealii::Point<dim, SizeType>     &q_point_loc,
                          const SizeType                         &element_volume,
                          Types::Index                            solve_block,
                          Types::Index index = Numbers::invalid_index) const = 0;

  /**
   * @brief User-implemented class for the RHS of postprocessed explicit equations.
   */
  virtual void
  compute_postprocess_explicit_rhs(VariableContainer<dim, degree, number> &variable_list,
                                   const SimulationTime               &simulation_time,
                                   const dealii::Point<dim, SizeType> &q_point_loc,
                                   const SizeType                     &element_volume,
                                   Types::Index solve_block) const = 0;

private:
  /**
   * @brief the timestep.
   */
  number delta_t;
};

PRISMS_PF_END_NAMESPACE