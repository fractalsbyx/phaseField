#ifndef SYSTEMCONTAINER_H
#define SYSTEMCONTAINER_H

#include <deal.II/base/exceptions.h>

#include "FieldContainer.h"
#include "ParaboloidSystem.h"

#include <core/userInputParameters.h>
#include <core/variableContainer.h>
#include <map>
#include <string>

template <int dim, int degree>
class SystemContainer
{
public:
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)

  struct PhaseData
  {
    FieldContainer<dim> omega;
    FieldContainer<dim> h;
  };

  struct CompData
  {
    FieldContainer<dim> mu;
    FieldContainer<dim> dmudt;
    scalarValue         M;
  };

  struct OPData
  {
    FieldContainer<dim>                        eta;
    FieldContainer<dim>                        detadt;
    std::map<std::string, FieldContainer<dim>> dhdeta;
  };

  const ParaboloidSystem         &isoSys;
  const userInputParameters<dim> &userInputs;

  std::map<std::string, PhaseData>            phase_data;
  std::map<std::string, CompData>             comp_data;
  std::vector<std::pair<std::string, OPData>> op_data;
  FieldContainer<dim>                         sum_sq_eta;

  /**
   * Constructor
   */
  SystemContainer(const ParaboloidSystem &sys, const userInputParameters<dim> &inputs)
    : isoSys(sys)
    , userInputs(inputs)
    , sum_sq_eta({}) // Zero initialize
  {
    // Initialize maps
    phase_data.clear();
    comp_data.clear();
    op_data.clear();
  }

  ~SystemContainer()
  {}

  void
  initialize_fields_explicit(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    for (const auto &comp_name : isoSys.comp_names)
      {
        comp_data[comp_name].mu.val  = variable_list.get_scalar_value(var_index);
        comp_data[comp_name].mu.grad = variable_list.get_scalar_gradient(var_index);
        var_index++;
      }
    for (const auto &phase_name : isoSys.order_params)
      {
        OPData op;
        op.eta.val  = variable_list.get_scalar_value(var_index);
        op.eta.grad = variable_list.get_scalar_gradient(var_index);
        op_data.push_back({phase_name, op});
        var_index++;
      }
  }

  void
  initialize_fields_postprocess(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    for (const auto &comp_name : isoSys.comp_names)
      {
        comp_data[comp_name].mu.val = variable_list.get_scalar_value(var_index);
        var_index++;
      }
    for (const auto &phase_name : isoSys.order_params)
      {
        OPData op;
        op.eta.val = variable_list.get_scalar_value(var_index);
        op_data.push_back({phase_name, op});
        var_index++;
      }
  }

  void
  calculate_omega_phase()
  {
    for (const auto &[phase_name, phase_info] : isoSys.phases)
      {
        PhaseData &phase = phase_data[phase_name];
        phase.omega.val  = phase_info.f_min;
        for (auto &[i_name, i_data] : comp_data)
          {
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              phase_info.comps.at(i_name);
            phase.omega +=
              -i_data.mu * i_data.mu / (2.0 * isoSys.Vm * isoSys.Vm * comp_info.k_well) -
              i_data.mu * comp_info.c_min / isoSys.Vm;
          }
      }
  }

  void
  calculate_sum_sq_eta()
  {
    sum_sq_eta.val = constV(0.);
    for (const auto &[phase_name, op] : op_data)
      {
        sum_sq_eta += op.eta * op.eta;
      }
  }

  void
  calculate_h()
  {
    for (auto &[phase_name, op] : op_data)
      {
        phase_data[phase_name].h += op.eta * op.eta;
      }
    for (auto &[phase_name, phase] : phase_data)
      {
        phase.h /= sum_sq_eta;
      }
  }

  void
  calculate_dhdeta()
  {
    for (auto &[alpha_name, op] : op_data)
      {
        for (auto &[beta_name, beta] : phase_data)
          {
            FieldContainer<dim> &dhdeta = op.dhdeta[beta_name];
            dhdeta.val                  = constV(0.);
            if (alpha_name == beta_name)
              {
                dhdeta += 2.0 * op.eta;
              }
            dhdeta -= 2.0 * op.eta * beta.h;
            dhdeta /= sum_sq_eta;
          }
      }
  }

  void
  calculate_detadt()
  {
    for (auto &[alpha_name, op] : op_data)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases.at(alpha_name);

        double m     = 6.00 * phase_info.sigma / isoSys.l_int;
        double kappa = 0.75 * phase_info.sigma * isoSys.l_int;
        double L     = 4.00 * phase_info.mu_int / isoSys.l_int / 3.00;

        // This is a variation, NOT a field
        FieldContainer<dim> interface_term;
        interface_term.val =
          m * (op.eta.val * op.eta.val * op.eta.val - op.eta.val +
               2. * 1.5 * op.eta.val * (sum_sq_eta.val - op.eta.val * op.eta.val));
        interface_term.grad = -kappa * op.eta.grad;

        // This is a variation, but has no vector term.
        scalarValue chemical_term = constV(0.);
        for (const auto &[beta_name, beta] : phase_data)
          {
            // NOTE: Only multiply values
            chemical_term += beta.omega.val * op.dhdeta.at(beta_name).val;
          }
        op.detadt = -L * (interface_term + chemical_term);
      }
  }

  void
  calculate_local_mobility()
  {
    for (auto &[comp_name, comp] : comp_data)
      {
        comp.M = constV(0.);
        for (auto &[phase_name, phase] : phase_data)
          {
            comp.M += isoSys.phases.at(phase_name).D * phase.h.val /
                      (isoSys.Vm * isoSys.Vm *
                       isoSys.phases.at(phase_name).comps.at(comp_name).k_well);
          }
      }
  }

  void
  calculate_dmudt()
  {
    for (auto &[comp_name, comp] : comp_data)
      {
        FieldContainer<dim> chi_AA;
        for (const auto &[phase_name, phase] : phase_data)
          {
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys.phases.at(phase_name).comps.at(comp_name);
            chi_AA += phase.h / comp_info.k_well / isoSys.Vm / isoSys.Vm;
          }
        comp.dmudt.val = constV(0.);
        // Flux term
        comp.dmudt.grad = -comp.M * -comp.mu.grad;

        // Partitioning term
        for (auto &[phase_name, op] : op_data)
          {
            FieldContainer<dim> drhodeta_sum;
            for (const auto &[beta_name, beta] : phase_data)
              {
                auto &comp_info = isoSys.phases.at(beta_name).comps.at(comp_name);
                drhodeta_sum +=
                  op.dhdeta.at(beta_name) *
                  (comp.mu / isoSys.Vm / comp_info.k_well + comp_info.c_min);
              }
            drhodeta_sum /= isoSys.Vm;
            comp.dmudt -= FieldContainer<dim>::field_x_variation(drhodeta_sum, op.detadt);
          }
        comp.dmudt = FieldContainer<dim>::field_x_variation(1.0 / chi_AA, comp.dmudt);
      }
  }

  void
  calculate_locals()
  {
    calculate_omega_phase();
    calculate_sum_sq_eta();
    calculate_h();
    calculate_dhdeta();
    calculate_local_mobility();
  }

  void
  submit_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                            &var_index)
  {
    for (auto &[comp_name, comp] : comp_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                comp.mu.val +
                                                  comp.dmudt.val * userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -comp.dmudt.grad * userInputs.dtValue);
        var_index++;
      }
    for (auto &[phase_name, op] : op_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                op.eta.val +
                                                  op.detadt.val * userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -op.detadt.grad * userInputs.dtValue);
        var_index++;
      }
  }

  void
  submit_pp_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
    uint                                                            &pp_index)
  {
    for (const auto &[comp_name, comp] : comp_data)
      {
        scalarValue c = constV(0.);
        for (const auto &[phase_name, phase] : phase_data)
          {
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys.phases.at(phase_name).comps.at(comp_name);
            c += phase.h.val * comp_info.c_min;
            c += phase.h.val * comp.mu.val / comp_info.k_well / isoSys.Vm;
          }
        pp_variable_list.set_scalar_value_term_RHS(pp_index, c);
        pp_index++;
      }
  }
};

#endif