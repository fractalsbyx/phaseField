#ifndef PARABOLOIDSYSTEM_H
#define PARABOLOIDSYSTEM_H

#include "json.hpp"

#include <core/variableAttributeLoader.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

class ParaboloidSystem
{
public:
  struct PhaseCompInfo
  {
    double c_min, k_well, x0;
    double _M, M;
  };

  struct Phase
  {
    double                               _mu_int, mu_int;
    double                               _D, D;
    double                               _sigma, sigma;
    double                               _f_min, f_min;
    std::map<std::string, PhaseCompInfo> comps;
  };

  std::map<std::string, Phase> phases;
  std::set<std::string>        comp_names;
  std::string                  solution_component;
  std::vector<std::string>     order_params;
  uint                         N;
  double                       length_scale, time_scale, energy_scale;
  double                       _Vm, Vm;
  double                       _l_int, l_int;

  ParaboloidSystem()
  {}

  ParaboloidSystem(const nlohmann::json &TCSystem)
  {
    from_json(TCSystem);
  }

  void
  from_json(const nlohmann::json &j)
  {
    // Parse Dimensions
    length_scale = j.at("dimensions").at("length_scale").get<double>();
    time_scale   = j.at("dimensions").at("time_scale").get<double>();
    energy_scale = j.at("dimensions").at("energy_density_scale").get<double>();

    // Parse Vm
    _Vm = j.at("Vm").get<double>();

    // Parse solution component
    solution_component = j.at("solution_component").get<std::string>();

    // Parse components
    for (const auto &comp_name : j.at("components"))
      {
        comp_names.insert(comp_name);
      }

    // Parse phases
    for (const auto &[phase_name, phase_data] : j.at("phases").items())
      {
        Phase phase;
        phase._mu_int = phase_data.at("mu_int").get<double>();
        phase._sigma  = phase_data.at("sigma").get<double>();

        for (const auto &[comp_name, comp_data] : phase_data.items())
          {
            if (comp_name != "mu_int" && comp_name != "sigma" && comp_name != "f_min" &&
                comp_name != "D")
              {
                PhaseCompInfo phaseCompInfo;
                phaseCompInfo.c_min    = comp_data.at("c_min").get<double>();
                phaseCompInfo.k_well   = comp_data.at("k_well").get<double>();
                phaseCompInfo.x0       = comp_data.at("x0").get<double>();
                phase.comps[comp_name] = phaseCompInfo;
              }
          }
        phases[phase_name] = phase;
      }
    N = phases.size();
    // Parse order parameters
    for (const auto &phase_name : j.at("order_parameters"))
      {
        order_params.push_back(phase_name);
      }
    nondimensionalize();
  }

  void
  nondimensionalize()
  {
    const double &l0 = length_scale;
    const double &t0 = time_scale;
    const double &E0 = energy_scale; // density
    Vm               = _Vm / (l0 * l0 * l0);
    l_int            = _l_int / (l0);
    for (auto &[phase_name, phase] : phases)
      {
        phase.mu_int = phase._mu_int * (E0 * t0);
        phase.D      = phase._D * (t0 / (l0 * l0)); // FIX?
        phase.sigma  = phase._sigma / (E0 * l0);
      }
  }

  void
  print_parameters()
  {
    // Set column width
    const int col_width = 15;

    // Print header
    std::cout << std::left << std::setw(col_width) << "Name" << std::setw(col_width)
              << "Dimensionless" << std::setw(col_width) << "Dimensional"
              << "\n";
    std::cout << std::string(45, '-') << "\n";

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
    for (const auto &[phase_name, phase] : phases)
      {
        std::cout << std::setw(col_width) << phase_name << "\n";
        std::cout << std::setw(col_width) << "mu_int:" << std::setw(col_width)
                  << phase.mu_int << std::setw(col_width) << phase._mu_int << "\n";
        std::cout << std::setw(col_width) << "D:" << std::setw(col_width) << phase.D
                  << std::setw(col_width) << phase._D << "\n";
        std::cout << std::setw(col_width) << "sigma:" << std::setw(col_width)
                  << phase.sigma << std::setw(col_width) << phase._sigma << "\n";

        for (const auto &[comp_name, comp] : phase.comps)
          {
            std::cout << std::setw(col_width) << comp_name << "\n";
            std::cout << std::setw(col_width) << "c_min:" << std::setw(col_width)
                      << comp.c_min << "\n";
            std::cout << std::setw(col_width) << "k_well:" << std::setw(col_width)
                      << comp.k_well << "\n";
            std::cout << std::setw(col_width) << "x0:" << std::setw(col_width) << comp.x0
                      << "\n";
          }
        std::cout << "\n";
      }
  }

  void
  load_variables(variableAttributeLoader *loader)
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
    std::vector<std::string>    op_names;
    std::vector<std::string>    grad_op_names;
    std::map<std::string, uint> phase_counter;
    for (const auto &phase_name : order_params)
      {
        if (phase_counter.find(phase_name) == phase_counter.end())
          {
            phase_counter[phase_name] = 0;
          }
        std::string var_name =
          phase_name + "_" + std::to_string(phase_counter[phase_name]);
        op_names.push_back(var_name);
        grad_op_names.push_back("grad(" + var_name + ")");
        phase_counter[phase_name]++;
        std::cout << var_name << " ";
      }
    std::cout << "\n";

    // Assign fields
    uint var_index = 0;
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

  // void
  // load_pp_variables(variableAttributeLoader *loader)
  //{}
};

#endif