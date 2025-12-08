#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/pde_problem.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>

#include <vector>

using namespace prisms;

int
main()
{
  // for my sanity
  const EvalFlags values               = EvalFlags::values;
  const EvalFlags gradients            = EvalFlags::gradients;
  const EvalFlags values_and_gradients = values | gradients;

  std::vector<FieldAttributes> fields;
  fields.reserve(7);

  fields.emplace_back("n0", Scalar, values_and_gradients);
  fields.emplace_back("n1", Scalar, values_and_gradients);
  fields.emplace_back("n2", Scalar, values_and_gradients);
  fields.emplace_back("c0", Scalar, values_and_gradients);
  fields.emplace_back("c1", Scalar, values_and_gradients);
  fields.emplace_back("phi", Scalar, values);
  fields.emplace_back("some_postprocess", Scalar, values);

  DependencySet op_vals_and_grads({
    {0, values_and_gradients},
    {1, values_and_gradients},
    {2, values_and_gradients}
  });
  DependencySet comp_vals({
    {3, values},
    {4, values}
  });

  SolveGroup order_parameters(1, ExplicitTimeDependent, {0, 1, 2});
  order_parameters.dependencies_rhs.insert(op_vals_and_grads.begin(),
                                           op_vals_and_grads.end());
  order_parameters.dependencies_rhs.insert(comp_vals.begin(), comp_vals.end());

  DependencySet       comp_vals_and_grads({
    {3, values_and_gradients},
    {4, values_and_gradients}
  });
  DependencySet       comp_new_and_old({
    {3, values_and_gradients, 1},
    {4, values_and_gradients, 1}
  });
  ChangeDependencySet comp_changes({
    {3, values_and_gradients},
    {4, values_and_gradients}
  });

  SolveGroup diffusion(2, ImplicitTimeDependent, {0, 1, 2});
  diffusion.dependencies_rhs.insert(comp_new_and_old.begin(), comp_new_and_old.end());
  diffusion.dependencies_rhs.insert(op_vals_and_grads.begin(), op_vals_and_grads.end());

  diffusion.dependencies_lhs.insert(comp_vals.begin(), comp_vals.end());
  diffusion.dependencies_lhs.insert(comp_vals_and_grads.begin(),
                                    comp_vals_and_grads.end());
  diffusion.dependencies_lhs.insert(op_vals_and_grads.begin(), op_vals_and_grads.end());
  diffusion.dependencies_change.insert(comp_changes.begin(), comp_changes.end());

  SolveGroup constant(3, Constant, {5});
  // Postprocess currently is not a PDEType
  // SolveGroup postprocess(4, Postprocess, {6});

  std::set<SolveGroup> equation_block({order_parameters, diffusion, constant});
  // ...

  PDEProblem<dim, degree, double> problem(fields,
                                          equation_block,
                                          user_inputs,
                                          pf_tools,
                                          pde_operator);
  problem.solve();
}
