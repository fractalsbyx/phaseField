// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{
#include "IsothermalSystem.h"
  std::cout << "Starting loadPostProcessorVariableAttributes...\n";
  // Open file ------------------------------------------------------------------------
  // Different scope from customPDE, must initialize identical system
  std::string sysFile = "system.json";
  // Create an input file stream
  std::ifstream inputFile(sysFile);
  // Check if the file was successfully opened
  if (!inputFile.is_open())
    {
      std::cerr << "Could not open the file: " << sysFile << std::endl;
      std::exit(1);
    }
  // Parse the JSON file into a JSON object
  nlohmann::json TCSystem;
  try
    {
      inputFile >> TCSystem;
    }
  catch (nlohmann::json::parse_error &e)
    {
      std::cerr << "JSON parse error: " << e.what() << std::endl;
      std::exit(1);
    }
  // Close the file -------------------------------------------------------------------
  inputFile.close();

  // Make System object
  IsothermalSystem Sys(TCSystem);

  // Load Variable Attributes
  Sys.load_pp_variables(this);

  std::cout << "Finished loadPostProcessorVariableAttributes\n";
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.h. It takes in "variable_list" and "q_point_loc" as inputs and outputs two
// terms in the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.h' for assigning the
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
  // --- Getting the values and derivatives of the model variables ---
  SystemContainer<dim, degree> sysFields(Sys, userInputs);
  // Solve
  // std::cout << "Initialize Fields Start...\n";
  sysFields.initialize_fields(variable_list);
  sysFields.submit_pp_fields(pp_variable_list);
}
