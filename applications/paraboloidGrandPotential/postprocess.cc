// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing
// variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but
// for the postprocessing expressions. It sets the attributes for each
// postprocessing expression, including its name, whether it is a vector or
// scalar (only scalars are supported at present), its dependencies on other
// variables and their derivatives, and whether to calculate an integral of the
// postprocessed quantity over the entire domain. Note: this function is not a
// member of customPDE.

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{
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
  uint pp_index = 0;
  isoSys.load_pp_variables(this, pp_index);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and
// 'nonExplicitEquationRHS' in equations.h. It takes in "variable_list" and
// "q_point_loc" as inputs and outputs two terms in the expression for the
// postprocessing variable -- one proportional to the test function and one
// proportional to the gradient of the test function. The index for each
// variable in this list corresponds to the index given at the top of this file
// (for submitting the terms) and the index in 'equations.h' for assigning the
// values/derivatives of the primary variables.

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  SystemContainer<dim, degree> ppsys(isoSys, userInputs);
  uint                         var_index = 0;
  ppsys.initialize_fields_postprocess(variable_list, var_index);
  ppsys.calculate_deltas();
  ppsys.calculate_sum_sq_eta();
  ppsys.calculate_h();
  uint pp_index = 0;
  ppsys.submit_pp_fields(pp_variable_list, pp_index);
}
