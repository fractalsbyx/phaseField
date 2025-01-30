// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for
// each function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void
customAttributeLoader::loadVariableAttributes()
{
// Load the model parameters
#include "ParaboloidSystem.h"

#include <fstream>
#include <iostream>
  nlohmann::json   model_parameters;
  ParaboloidSystem isoSys;
  std::ifstream    system_file;
  system_file.open("system.json");
  system_file >> model_parameters;
  system_file.close();
  isoSys.from_json(model_parameters);
  isoSys.load_variables(this);
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a
// list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one
// proportional to the test function and one proportional to the gradient of the
// test function. The index for each variable in this list corresponds to the
// index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  SystemContainer<dim, degree> sys(isoSys, userInputs);
  uint                         var_index = 0;
  sys.initialize_fields_explicit(variable_list, var_index);
  sys.calculate_locals();
  sys.calculate_detadt();
  sys.calculate_dmudt();
  var_index = 0;
  sys.submit_fields(variable_list, var_index);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are
// not explicit time-dependent equations. It takes "variable_list" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is
// given by "q_point_loc". The function outputs two terms to variable_list --
// one proportional to the test function and one proportional to the gradient of
// the test function. The index for each variable in this list corresponds to
// the index given at the top of this file.

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. The
// (x,y,z) location of that quadrature point is given by "q_point_loc". The
// function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function -- for the
// left-hand-side of the equation. The index for each variable in this list
// corresponds to the index given at the top of this file. If there are multiple
// elliptic equations, conditional statements should be sed to ensure that the
// correct residual is being submitted. The index of the field being solved can
// be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}

#include "SystemContainer.h"
// template class SystemContainer<1, 1>;
// template class SystemContainer<1, 2>;
// template class SystemContainer<1, 3>;
// template class SystemContainer<1, 4>;
// template class SystemContainer<1, 5>;
// template class SystemContainer<1, 6>;

template class SystemContainer<2, 1>;
template class SystemContainer<2, 2>;
template class SystemContainer<2, 3>;
template class SystemContainer<2, 4>;
template class SystemContainer<2, 5>;
template class SystemContainer<2, 6>;

template class SystemContainer<3, 1>;
template class SystemContainer<3, 2>;
template class SystemContainer<3, 3>;
template class SystemContainer<3, 4>;
template class SystemContainer<3, 5>;
template class SystemContainer<3, 6>;

template class customPDE<2, 1>;
template class customPDE<2, 2>;
template class customPDE<2, 3>;
template class customPDE<2, 4>;
template class customPDE<2, 5>;
template class customPDE<2, 6>;

template class customPDE<3, 1>;
template class customPDE<3, 2>;
template class customPDE<3, 3>;
template class customPDE<3, 4>;
template class customPDE<3, 5>;
template class customPDE<3, 6>;