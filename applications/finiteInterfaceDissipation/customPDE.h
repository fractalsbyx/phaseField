#include "SystemContainer.h"
#include "matrixFreePDE.h"

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {
    std::cout << "Starting customPDE initializer...\n";
    // Read system json file to TCSystem
    parseSystem();
    Sys.print_parameters();
    print_initial_energies();
    std::cout << "customPDE initialized\n";
    for (uint i = 0; i < this->var_attributes.attributes.size(); i++)
      {
        const auto &_set = this->var_attributes.attributes.at(i).dependencies_RHS;
        for (const auto &_var : _set)
          {
            std::cout << std::to_string(i) << ": " << _var << " "
                      << this->var_attributes.attributes.at(i).is_nonlinear << "\n";
          }
      }
  };

  void
  print_initial_energies()
  {
    SystemContainer<dim, degree> sys_for_print(Sys, userInputs);
    for (const auto &[phase_name, phase] : Sys.phases)
      {
        for (const auto &[comp_name, comp_info] : phase.comps)
          {
            sys_for_print.phase_fields[phase_name]->comp_data[comp_name].x_data.val =
              constV(comp_info.x0);
          }
      }
    sys_for_print.calculate_locals();
    std::cout << "Initial bulk free energies:\n";
    for (const auto &[phase_name, phase] : Sys.phases)
      {
        std::cout << phase_name << ":\n"
                  << "Free Energy:\t"
                  << sys_for_print.phase_fields[phase_name]->phase_free_energy[0] << "\n";
        for (const auto &[comp_name, comp_info] : phase.comps)
          {
            std::cout
              << "dfdx_" << comp_name << ":\t"
              << sys_for_print.phase_fields[phase_name]->comp_data[comp_name].dfdx.val[0]
              << "\n";
          }
        std::cout << "\n";
      }
  }

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include "typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.cc)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.cc)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the LHS of the governing equations (in equations.cc)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

// Function to set postprocessing expressions (in postprocess.h)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif

// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================

  void
  parseSystem()
  {
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
    // Close the file
    inputFile.close();

    Sys = IsothermalSystem(TCSystem);
  }

  inline double
  interface(const double x)
  {
    return 0.5 * (1.0 + std::sin(pi * std::max(-0.5, std::min(0.5, x / Sys.eta))));
  }

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================

  IsothermalSystem Sys;
  // File for system name
  std::string sysFile = "system.json"; // userInputs.get_model_constant_string("sysFile");

  double T;
  double sigma;
  double l_gb;

  // prm constants
  double r0 = userInputs.get_model_constant_double("r0");

  // Zero vector
  scalargradType ZERO;
  // pi
  const double pi = 3.141592653589793238;

  // ================================================================
};
