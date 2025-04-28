#ifndef PARABOLOIDSYSTEM_H
#define PARABOLOIDSYSTEM_H

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "json.hpp"

#include <core/variableAttributeLoader.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

template <typename number>
using boost_vector = boost::numeric::ublas::vector<number>;
template <typename number>
using boost_matrix = boost::numeric::ublas::symmetric_matrix<number>;

/**
 * @brief Class containing the thermodynamic and kinetic parameters needed for the grand
 * potential based model using second-order polynomials
 */
class ParaboloidSystem
{
public:
  /**
   * @brief Parameters for each phase
   */
  struct Phase
  {
    std::string name;
    double      _mu_int, mu_int;
    double      _D, D;
    double      _sigma, sigma;

    boost_matrix<double>       _A_well, A_well, suscept;
    boost_vector<double>       _B_well, B_well, c_ref, c0, mu0;
    double                     _D_well, D_well;
    dealii::FullMatrix<double> helper, helper_inv;
  };

  /**
   * @brief Phase parameters
   */
  std::vector<Phase> phases;
  /**
   * @brief Component names at each index
   */
  std::vector<std::string> comp_names;
  /**
   * @brief Number of non-solution components. (Binary->1, Ternary->2...)
   */
  unsigned int num_comps;
  /**
   * @brief Phase names at each index
   */
  std::vector<std::string> phase_names;
  /**
   * @brief Name of the component that is not an independent variable
   */
  std::string solution_component;
  /**
   * @brief The index of the phase of each order parameter
   */
  std::vector<uint> order_params;
  /**
   * @brief Scale parameters for non-dimensionalization
   * - length_scale: length scale
   * - time_scale: time scale
   * - energy_scale: energy density scale
   */
  double length_scale, time_scale, energy_scale;
  /**
   * @brief The atomic or molar volume assumed to be uniform
   */
  double _Vm, Vm;
  /**
   * @brief The interface width
   */
  double _l_int, l_int;
  /**
   * @brief If true, the energy density is converted to volumetric energy density from
   * molar/atomic energy
   */
  bool volumetrize;

  /**
   * @brief Constructor
   */
  ParaboloidSystem()
  {}

  /**
   * @brief JSON Constructor
   * @param TCSystem JSON object containing the parameters
   */
  ParaboloidSystem(const nlohmann::json &TCSystem)
  {
    from_json(TCSystem);
  }

  /**
   * @brief Load the parameters from a JSON object
   * @param j JSON object containing the parameters
   */
  void
  from_json(const nlohmann::json &j)
  {
    // Parse Dimensions
    length_scale = j.at("dimensions").at("length_scale").get<double>();
    time_scale   = j.at("dimensions").at("time_scale").get<double>();
    energy_scale = j.at("dimensions").at("energy_density_scale").get<double>();

    // Parse Vm
    _Vm    = j.at("Vm").get<double>();
    _l_int = j.at("l_int").get<double>();

    // Check if volumetrization needed
    volumetrize = j.at("convert_fractional_to_volumetric_energy").get<bool>();

    // Parse solution component
    solution_component = j.at("solution_component").get<std::string>();

    // Parse components
    comp_names.clear();
    for (const auto &comp_name : j.at("components"))
      {
        comp_names.push_back(comp_name);
      }
    num_comps = comp_names.size();

    // Parse phases
    phases.clear();
    phase_names.clear();
    for (const auto &[phase_name, phase_data] : j.at("phases").items())
      {
        Phase phase;
        phase.name    = phase_name;
        phase._mu_int = phase_data.at("mu_int").get<double>();
        phase._sigma  = phase_data.at("sigma").get<double>();
        phase._D      = phase_data.at("D").get<double>();

        phase._A_well = boost_matrix<double>(num_comps, num_comps);
        phase._B_well = boost_vector<double>(num_comps);
        phase._D_well = phase_data.at("D_well").get<double>();

        // Parse components
        for (uint comp_index = 0; comp_index < num_comps; comp_index++)
          {
            const std::string &comp_name = comp_names[comp_index];
            phase.c_ref[comp_index] = phase_data.at("c_ref").at(comp_name).get<double>();
            phase.B_well[comp_index] =
              phase_data.at("B_well").at(comp_name).get<double>();
            phase.c0[comp_index] = phase_data.at("c0").at(comp_name).get<double>();

            auto A_well_row = phase_data.at(comp_name).at("A_well");
            for (uint col_comp_index = 0; col_comp_index < num_comps; col_comp_index++)
              {
                const std::string &col_comp_name = comp_names[col_comp_index];
                if (!A_well_row.contains(col_comp_name))
                  {
                    continue;
                  }
                phase._A_well(comp_index, col_comp_index) =
                  A_well_row.at(col_comp_name).get<double>();
              }
          }
        phases.push_back(phase);
        phase_names.push_back(phase_name);
      }

    // Parse order parameters
    order_params.clear();
    for (const auto &phase_name : j.at("order_parameters"))
      {
        uint phase_index = std::find(phase_names.begin(), phase_names.end(), phase_name) -
                           phase_names.begin();
        order_params.push_back(phase_index);
      }

    // Convert to volumetric energy if necessary
    if (volumetrize)
      {
        for (Phase &phase : phases)
          {
            phase._A_well /= _Vm;
            phase._B_well /= _Vm;
            phase._D_well /= _Vm;
          }
      }
    // Non-dimensionalize
    nondimensionalize();
  }

  /**
   * @brief Non-dimensionalize the parameters using the provided unit scales.
   */
  void
  nondimensionalize()
  {
    const double &l0 = length_scale;
    const double &t0 = time_scale;
    const double &E0 = energy_scale; // energy density
    Vm               = _Vm / (l0 * l0 * l0);
    l_int            = _l_int / (l0);
    for (Phase &phase : phases)
      {
        phase.mu_int = phase._mu_int / (l0 / (E0 * t0));
        phase.sigma  = phase._sigma / (E0 * l0);
        phase.D      = phase._D / ((l0 * l0) / t0);

        phase.A_well = phase._A_well;
        phase.A_well /= (E0);
        phase.B_well = phase._B_well;
        phase.B_well /= (E0);
        phase.D_well = phase._D_well;
        phase.D_well /= (E0);

        phase.mu0 = prod(phase.A_well, phase.c0) + phase.B_well;

        phase.helper = dealii::FullMatrix<double>(num_comps);
        for (size_t i = 0; i < num_comps; ++i)
          {
            for (size_t j = 0; j < num_comps; ++j)
              {
                phase.helper(i, j) = phase.A_well(i, j);
              }
          }
        phase.helper_inv.invert(phase.helper);
        phase.suscept = boost_matrix<double>(num_comps, num_comps);
        for (size_t i = 0; i < num_comps; ++i)
          {
            for (size_t j = 0; j < num_comps; ++j)
              {
                phase.suscept(i, j) = phase.helper_inv(i, j);
              }
          }
      }
  }

  /**
   * @brief Print the parameters to the console
   */
  void
  print_parameters()
  {
    // Set column width
    const int col_width  = 15;
    const int line_width = 45;

    // Print header
    std::cout << std::left << std::setw(col_width) << "Name" << std::setw(col_width)
              << "Dimensionless" << std::setw(col_width) << "Dimensional"
              << "\n";
    std::cout << std::string(line_width, '-') << "\n";

    // Print Vm and l_int
    std::cout << std::setw(col_width) << "Vm:" << std::setw(col_width) << Vm
              << std::setw(col_width) << _Vm << "\n";
    std::cout << std::setw(col_width) << "l_int:" << std::setw(col_width) << l_int
              << std::setw(col_width) << _l_int << "\n";

    // Print component information
    for (const auto &comp_name : comp_names)
      {
        std::cout << std::setw(col_width) << comp_name << "\n";
      }
    std::cout << "\n";

    // Print phase information
    for (const Phase &phase : phases)
      {
        std::cout << std::setw(col_width) << phase.name << "\n";
        std::cout << std::setw(col_width) << "mu_int:" << std::setw(col_width)
                  << phase.mu_int << std::setw(col_width) << phase._mu_int << "\n";
        std::cout << std::setw(col_width) << "D:" << std::setw(col_width) << phase.D
                  << std::setw(col_width) << phase._D << "\n";
        std::cout << std::setw(col_width) << "sigma:" << std::setw(col_width)
                  << phase.sigma << std::setw(col_width) << phase._sigma << "\n";

        // Print A matrix
        std::cout << std::setw(col_width) << "A_well:\n";
        std::cout << std::setw(col_width) << "..";
        for (uint col_index = 0; col_index < comp_names.size(); col_index++)
          {
            std::cout << std::setw(col_width) << comp_names[col_index];
          }
        for (uint row_index = 0; row_index < comp_names.size(); row_index++)
          {
            std::cout << std::setw(col_width) << comp_names[row_index];
            for (uint col_index = 0; col_index < comp_names.size(); col_index++)
              {
                std::cout << std::setw(col_width) << phase.A_well(row_index, col_index);
              }
            std::cout << "\n";
          }
        std::cout << "\n";

        // Print B vector
        std::cout << std::setw(col_width) << "B_well:\n";
        // std::cout << std::setw(col_width) << "..";
        for (uint comp_index = 0; comp_index < comp_names.size(); comp_index++)
          {
            std::cout << std::setw(col_width) << comp_names[comp_index];
          }
        std::cout << "\n";
        for (uint comp_index = 0; comp_index < comp_names.size(); comp_index++)
          {
            std::cout << std::setw(col_width) << phase.B_well[comp_index];
          }
        std::cout << "\n";

        // Print D_well
        std::cout << std::setw(col_width) << "D_well:" << std::setw(col_width)
                  << phase.D_well << std::setw(col_width) << phase._D_well << "\n";
        std::cout << "\n";
      }
  }

  /**
   * @brief Declare the fields needed for the PDE in PRISMS-PF
   * @param loader Pointer to the attribute loader (equations.cc)
   */
  void
  load_variables(customAttributeLoader *loader, uint &var_index)
  {
    // Get names for mu fields
    std::vector<std::string> mu_names;
    std::vector<std::string> grad_mu_names;
    std::cout << "Comp names: ";
    for (const auto &comp_name : comp_names)
      {
        std::string var_name = "mu_" + comp_name;
        mu_names.push_back(var_name);
        grad_mu_names.push_back("grad(" + var_name + ")");
        std::cout << var_name << " ";
      }
    std::cout << "\n"
              << "Order parameter names: ";
    // Get names for order parameter fields
    std::vector<std::string> op_names;
    std::vector<std::string> grad_op_names;
    std::map<uint, uint>     phase_counter;
    for (const auto &phase_index : order_params)
      {
        if (phase_counter.find(phase_index) == phase_counter.end())
          {
            phase_counter[phase_index] = 0;
          }
        std::string phase_name = phase_names[phase_index];
        std::string var_name =
          phase_name + "_" + std::to_string(phase_counter[phase_index]);
        op_names.push_back(var_name);
        grad_op_names.push_back("grad(" + var_name + ")");
        phase_counter[phase_index]++;
        std::cout << var_name << " ";
      }
    std::cout << "\n";

    // Assign fields
    for (const auto &mu_name : mu_names)
      {
        loader->set_variable_name(var_index, mu_name);
        loader->set_variable_type(var_index, SCALAR);
        loader->set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);
        // loader->set_allowed_to_nucleate			(var_index, (phase_index>0));
        loader->set_need_value_nucleation(var_index, true);
        loader->insert_dependencies_value_term_RHS(var_index, mu_names);
        loader->insert_dependencies_value_term_RHS(var_index, grad_mu_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, mu_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, grad_mu_names);
        var_index++;
      }
    for (const auto &op_name : op_names)
      {
        loader->set_variable_name(var_index, op_name);
        loader->set_variable_type(var_index, SCALAR);
        loader->set_variable_equation_type(var_index, EXPLICIT_TIME_DEPENDENT);
        loader->set_need_value_nucleation(var_index, true);
        loader->insert_dependencies_value_term_RHS(var_index, mu_names);
        loader->insert_dependencies_value_term_RHS(var_index, grad_mu_names);
        loader->insert_dependencies_value_term_RHS(var_index, op_names);
        loader->insert_dependencies_value_term_RHS(var_index, grad_op_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, mu_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, grad_mu_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, op_names);
        loader->insert_dependencies_gradient_term_RHS(var_index, grad_op_names);
        var_index++;
      }
  }

  /**
   * @brief Declare the post-processed fields in PRISMS-PF
   * @param loader Pointer to the attribute loader (postprocess.cc)
   */
  void
  load_pp_variables(customAttributeLoader *loader, uint &pp_index)
  {
    // Get names for comp fields
    std::vector<std::string> c_names;
    std::vector<std::string> mu_names;
    for (const auto &comp_name : comp_names)
      {
        c_names.push_back("c_" + comp_name);
        mu_names.push_back("mu_" + comp_name);
      }
    // Get names for order parameter fields
    std::vector<std::string> op_names;
    std::map<uint, uint>     phase_counter;
    for (const auto &phase_index : order_params)
      {
        if (phase_counter.find(phase_index) == phase_counter.end())
          {
            phase_counter[phase_index] = 0;
          }
        std::string phase_name = phase_names[phase_index];
        std::string var_name =
          phase_name + "_" + std::to_string(phase_counter[phase_index]);
        op_names.push_back(var_name);
        phase_counter[phase_index]++;
      }
    std::cout << "\n";

    // Assign fields
    for (const auto &c_name : c_names)
      {
        loader->set_variable_name(pp_index, c_name);
        loader->set_variable_type(pp_index, SCALAR);
        loader->set_variable_equation_type(pp_index, EXPLICIT_TIME_DEPENDENT);
        loader->insert_dependencies_value_term_RHS(pp_index, mu_names);
        loader->insert_dependencies_value_term_RHS(pp_index, op_names);
        pp_index++;
      }
  }
};

#endif