// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================

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
  uint var_index = 0;
  isoSys.print_parameters();
  isoSys.load_variables(this, var_index);
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================

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

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}