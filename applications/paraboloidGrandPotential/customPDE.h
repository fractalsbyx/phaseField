#include "ParaboloidSystem.h"
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
    print_interface_properties();
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
  /**
   * @brief Function to set the initial conditions with the proper tanh profile given a
   * level-set function
   */
  [[nodiscard]] double
  interface(double x) const
  {
    return 0.5 * (1.0 + std::tanh(2.0 * x / isoSys.l_int));
  }

  /**
   * @brief Prints the grand potential densities (nondimensionalized) of each phase at
   * its initial conditions
   */
  void
  print_initial_energies()
  {
    std::cout << "Initial omega free energies:\n";
    for (uint phase_index = 0; phase_index < isoSys.phases.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase = isoSys.phases.at(phase_index);
        SystemContainer<dim, degree>   sys_for_print(isoSys, userInputs);
        sys_for_print.op_data.push_back({
          phase_index,
          {
            {constV(1.), {}}, // eta
            {constV(0.), {}}, // detadt
            {}                // dhdeta
          }
        });

        boost_vector<double> mu0 = prod(phase.A_well, phase.c0);
        for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
          {
            sys_for_print.mu[comp_index].val = mu0[comp_index];
          }
        sys_for_print.calculate_deltas();
        sys_for_print.calculate_omega_phase();

        std::cout << phase.name << ":\n"
                  << "Omega:\t" << sys_for_print.phase_data[phase_index].omega.val[0]
                  << "\n";
        initial_omega[phase.name] = sys_for_print.phase_data[phase_index].omega.val[0];
        for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
          {
            std::cout << "mu_" << isoSys.comp_names.at(comp_index) << ":\t"
                      << mu0[comp_index] << "\n";
          }
        std::cout << "\n";
      }
  }

  /**
   * @brief Function to estimate the numerical stability of the system based on the
   * parameters provided and give a recommendation for the time step
   */
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
    std::string      limiting_phase;
    for (const auto &phase : isoSys.phases)
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
            limiting_phase      = phase.name;
          }
        if (order_parameter_gradient_factor > max_gradient_factor)
          {
            max_gradient_factor = order_parameter_gradient_factor;
            stability_limiter   = "order parameter evolution";
            limiting_phase      = phase.name;
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

  /**
   * @brief Function to print the properties of the interface between two phases at
   * initial conditions
   */
  void
  print_interface_properties()
  {
    for (const auto &alpha : isoSys.phases)
      {
        for (const auto &beta : isoSys.phases)
          {
            std::cout << "Properties of the interface between " << alpha.name << " and "
                      << beta.name << ":\n";
            double delta_g = initial_omega[alpha.name] - initial_omega[beta.name];
            double sigma   = 0.5 * (alpha.sigma + beta.sigma);
            double D       = 0.5 * (alpha.D * beta.D) / (alpha.D + beta.D);
            double mu_int =
              0.5 * (alpha.mu_int * beta.mu_int) / (alpha.mu_int + beta.mu_int);
            std::cout << "The value of the dimensionless number delta_g/(sigma/l_int) is "
                      << delta_g * isoSys.l_int / sigma << "\n";
            std::cout
              << "The value of the dimensionless number delta_g*mu_int*l_int/D is "
              << delta_g * mu_int * isoSys.l_int / D << "\n";
            std::cout << "The value of the dimensionless number sigma*mu_int/D is "
                      << sigma * mu_int / D << "\n\n";
          }
      }
  }

  // ================================================================
  // Members specific to this subclass
  // ================================================================
  /**
   * @brief JSON containing the model parameters
   */
  nlohmann::json model_parameters;
  /**
   * @brief Object containing the thermodynamic and kinetic parameters
   */
  ParaboloidSystem isoSys;
  /**
   * @brief Fraction of the theoretical maximum time step to use
   */
  double timestep_alpha = userInputs.get_model_constant_double("timestep_alpha");
  /**
   * @brief Map of the initial grand potential densities for each phase (used for
   * printing)
   */
  std::map<std::string, double> initial_omega;

  double r0 = userInputs.get_model_constant_double("r0");
  // ================================================================
};
