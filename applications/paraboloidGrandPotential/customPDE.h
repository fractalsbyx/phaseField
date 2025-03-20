#include "SystemContainer.h"

#include <core/matrixFreePDE.h>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>

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
    // Load the model parameters
    std::ifstream ifs("system.json");
    ifs >> model_parameters;
    ifs.close();
    isoSys.from_json(model_parameters);
    isoSys.print_parameters();
    print_initial_energies();
    estimate_stability();
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
#include <core/typeDefs.h>

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the LHS of the governing equations (in equations.h)
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
  [[nodiscard]] double
  interface(double x) const
  {
    return 0.5 * (1.0 + std::tanh(2.0 * x / isoSys.l_int));
  }

  void
  print_initial_energies()
  {
    std::cout << "Initial omega free energies:\n";
    for (const auto &[phase_name, phase] : isoSys.phases)
      {
        SystemContainer<dim, degree> sys_for_print(isoSys, userInputs);
        sys_for_print.op_data.push_back({
          phase_name,
          {
            {constV(1.), {}}, // eta
            {constV(0.), {}}, // detadt
            {}                // dhdeta
          }
        });

        for (const auto &[comp_name, comp_info] : isoSys.phases.at(phase_name).comps)
          {
            double mu0 = isoSys.Vm * comp_info.k_well * (comp_info.x0 - comp_info.c_min);
            sys_for_print.comp_data[comp_name].mu.val = constV(mu0);
          }

        sys_for_print.calculate_locals();

        std::cout << phase_name << ":\n"
                  << "Omega:\t" << sys_for_print.phase_data[phase_name].omega.val[0]
                  << "\n";
        for (const auto &[comp_name, comp_info] : phase.comps)
          {
            std::cout << "mu_" << comp_name << ":\t"
                      << sys_for_print.comp_data[comp_name].mu.val[0] << "\n";
          }
        std::cout << "\n";
      }
  }

  void
  estimate_stability()
  {
    double min_dx = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < dim; i++)
      {
        min_dx = std::min(min_dx,
                          double(userInputs.subdivisions[i]) * userInputs.domain_size[i] /
                            std::pow(2.0, userInputs.refine_factor));
      }
    constexpr double theoretical_max_gradient_factor = 0.25;
    double           max_gradient_factor             = 0.0;
    std::string      stability_limiter               = "diffusion";
    std::string      limiting_phase                  = "";
    for (const auto &[phase_name, phase] : isoSys.phases)
      {
        const double gradient_prefactor = userInputs.dtValue *
                                          (userInputs.degree * userInputs.degree) /
                                          (min_dx * min_dx);

        const double diffusion_gradient_factor = phase.D * gradient_prefactor;
        // mu*sigma = L*kappa
        const double order_parameter_gradient_factor =
          (phase.mu_int * phase.sigma) * gradient_prefactor;

        if (diffusion_gradient_factor > max_gradient_factor)
          {
            max_gradient_factor = diffusion_gradient_factor;
            stability_limiter   = "diffusion";
            limiting_phase      = phase_name;
          }
        if (order_parameter_gradient_factor > max_gradient_factor)
          {
            max_gradient_factor = order_parameter_gradient_factor;
            stability_limiter   = "order parameter evolution";
            limiting_phase      = phase_name;
          }
      }
    std::cout << "Then numerical stability for this set of parameters is limited by "
              << stability_limiter << " in phase " << limiting_phase << ".\n";
    std::cout << "The reccommended maximum time step (using a design factor of "
              << timestep_alpha
              << ") is " /* << std::scientific << std::setprecision(2) */
              << timestep_alpha * userInputs.dtValue * theoretical_max_gradient_factor /
                   max_gradient_factor
              << ".\n\n" /* << std::defaultfloat << std::setprecision(6) */;
  }

  // ================================================================
  // Model constants specific to this subclass
  // ================================================================
  nlohmann::json   model_parameters;
  ParaboloidSystem isoSys;
  double           r0   = userInputs.get_model_constant_double("r0");
  double timestep_alpha = userInputs.get_model_constant_double("timestep_alpha");

  // ================================================================
};
