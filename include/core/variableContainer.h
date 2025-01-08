// This class permits the access of a subset of indexed fields and gives an
// error if any non-allowed fields are requested
#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "core/solveTypeEnums.h"
#include "core/variableAttributes.h"
#include "varTypeEnums.h"

#include <core/model_variables.h>

template <int dim, int degree, typename T>
class variableContainer
{
public:
#include <core/typeDefs.h>
  // Constructors

  // Standard contructor, used for most situations
  variableContainer(const dealii::MatrixFree<dim, double> &data,
                    AttributesList                         _src_attributes,
                    AttributesList                         _dst_attributes,
                    solveType                              _solve_type);

  // Methods to get the value/grad/hess in the residual method (this is how the
  // user gets these values in equations.h)
  T
  get_scalar_value(unsigned int global_variable_index) const;
  dealii::Tensor<1, dim, T>
  get_scalar_gradient(unsigned int global_variable_index) const;
  dealii::Tensor<2, dim, T>
  get_scalar_hessian(unsigned int global_variable_index) const;
  dealii::Tensor<1, dim, T>
  get_vector_value(unsigned int global_variable_index) const;
  dealii::Tensor<2, dim, T>
  get_vector_gradient(unsigned int global_variable_index) const;
  dealii::Tensor<3, dim, T>
  get_vector_hessian(unsigned int global_variable_index) const;

  T
  get_change_in_scalar_value(unsigned int global_variable_index) const;
  dealii::Tensor<1, dim, T>
  get_change_in_scalar_gradient(unsigned int global_variable_index) const;
  dealii::Tensor<2, dim, T>
  get_change_in_scalar_hessian(unsigned int global_variable_index) const;
  dealii::Tensor<1, dim, T>
  get_change_in_vector_value(unsigned int global_variable_index) const;
  dealii::Tensor<2, dim, T>
  get_change_in_vector_gradient(unsigned int global_variable_index) const;
  dealii::Tensor<3, dim, T>
  get_change_in_vector_hessian(unsigned int global_variable_index) const;

  // Methods to set the value residual and the gradient residual (this is how
  // the user sets these values in equations.h)
  void
  set_scalar_value_term_RHS(unsigned int global_variable_index, T val);
  void
  set_scalar_gradient_term_RHS(unsigned int              global_variable_index,
                               dealii::Tensor<1, dim, T> grad);
  void
  set_vector_value_term_RHS(unsigned int              global_variable_index,
                            dealii::Tensor<1, dim, T> val);
  void
  set_vector_gradient_term_RHS(unsigned int              global_variable_index,
                               dealii::Tensor<2, dim, T> grad);

  void
  set_scalar_value_term_LHS(unsigned int global_variable_index, T val);
  void
  set_scalar_gradient_term_LHS(unsigned int              global_variable_index,
                               dealii::Tensor<1, dim, T> grad);
  void
  set_vector_value_term_LHS(unsigned int              global_variable_index,
                            dealii::Tensor<1, dim, T> val);
  void
  set_vector_gradient_term_LHS(unsigned int              global_variable_index,
                               dealii::Tensor<2, dim, T> grad);

  // Initialize, read DOFs, and set evaulation flags for each variable
  void
  reinit_and_eval(const std::vector<vectorType *> &src, unsigned int cell);
  void
  reinit_and_eval_change_in_solution(const vectorType &src,
                                     unsigned int      cell,
                                     const uint       &var_index);

  /**
   * \brief Initialize the FEEvals for the destination variables
   */
  void
  reinit_dst(unsigned int cell);

  // Integrate the residuals and distribute from local to global
  void
  integrate_and_distribute(std::vector<vectorType *> &dst);
  void
  integrate_and_distribute_change_in_solution_LHS(vectorType &dst, const uint &var_index);

  // The quadrature point index, a method to get the number of quadrature points
  // per cell, and a method to get the xyz coordinates for the quadrature point
  unsigned int q_point = 0;

  [[nodiscard]] unsigned int
  get_num_q_points() const;

  [[nodiscard]] dealii::Point<dim, T>
  get_q_point_location() const;

private:
  // Object containing some information about each variable (indices, whether
  // the val/grad/hess is needed, etc)
  const AttributesList src_attributes;
  const AttributesList dst_attributes;
  const solveType      solve_type;

  // Vectors of the actual FEEvaluation objects for each active variable, split
  // into scalar variables and vector variables for type reasons
  using scalar_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, 1, double>;
  using vector_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, dim, double>;

  std::unordered_map<unsigned int, scalar_FEEval> src_scalar_vars_map;
  std::unordered_map<unsigned int, vector_FEEval> src_vector_vars_map;
  std::unordered_map<unsigned int, scalar_FEEval> dst_scalar_vars_map;
  std::unordered_map<unsigned int, vector_FEEval> dst_vector_vars_map;

  std::unordered_map<unsigned int, scalar_FEEval> scalar_change_in_vars_map;
  std::unordered_map<unsigned int, vector_FEEval> vector_change_in_vars_map;

  void
  AssertValid(const uint      &var_index,
              const fieldType &field_type,
              const EvalFlags &eval_flag,
              bool             is_change) const;
};

#endif
