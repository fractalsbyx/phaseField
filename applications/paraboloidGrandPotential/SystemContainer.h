#ifndef SYSTEMCONTAINER_H
#define SYSTEMCONTAINER_H

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
    FieldContainer<dim> phi;
    FieldContainer<dim> omega;
    FieldContainer<dim> h;
  };

  struct CompData
  {
    FieldContainer<dim> mu;
    FieldContainer<dim> dmudt;
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
  FieldContainer<dim>                         D;

  /**
   * Constructor
   */
  SystemContainer(const ParaboloidSystem         &_isoSys,
                  const userInputParameters<dim> &_userInputs)
    : isoSys(_isoSys)
    , userInputs(_userInputs)
  {}

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
  calculate_omega_phase()
  {
    for (const auto &[phase_name, phase_info] : isoSys.phases)
      {
        PhaseData &phase = phase_data[phase_name];
        phase.omega.val  = phase_info.f_min;
        for (auto &[i_name, i_data] : comp_data)
          {
            phase.omega += -i_data.mu * i_data.mu / (2.0 * isoSys.Vm * isoSys.Vm) //
                           - i_data.mu * phase_info.comps.at(i_name).c_min / isoSys.Vm;
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
                dhdeta += op.eta * 2.0;
              }
            dhdeta -= op.eta * beta.h * 2.0;
            dhdeta /= sum_sq_eta;
          }
      }
  }

  void
  calculate_detadt()
  {
    for (auto &[phase_name, op] : op_data)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases.at(phase_name);

        double m     = 6.00 * phase_info.sigma / isoSys.l_int;
        double kappa = 0.75 * phase_info.sigma * isoSys.l_int;
        double L     = 4.00 * phase_info.mu_int / isoSys.l_int / 3.00;

        // This is a variation, NOT a field
        FieldContainer<dim> interface_term;
        interface_term.val =
          m * (op.eta.val * op.eta.val * op.eta.val - op.eta.val +
               2. * 1.5 * op.eta.val * (sum_sq_eta.val - op.eta.val * op.eta.val));
        interface_term.grad = kappa * op.eta.grad;

        // This is a variation, but has no vector term.
        scalarValue chemical_term;
        for (auto &[phase_name, phase] : phase_data)
          {
            // NOTE: Only multiply values
            chemical_term += phase.omega.val * op.dhdeta.at(phase_name).val;
          }
        op.detadt = (interface_term + chemical_term) * (-L);
      }
  }

  void
  calculate_local_diffusivity()
  {
    D.val = constV(0.);
    for (auto &[phase_name, phase] : phase_data)
      {
        D += phase.h * isoSys.phases.at(phase_name).D;
      }
  }

  void
  calculate_dmudt()
  {
    for (auto &[comp_name, comp] : comp_data)
      {
        comp.dmudt.val = constV(0.);
        // Flux term
        comp.dmudt.grad = -D.val * comp.mu.grad; // CHECK SIGN

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
            comp.dmudt -=
              FieldContainer<dim>::field_x_variation(drhodeta_sum, op.detadt) *
              isoSys.Vm * isoSys.Vm *
              isoSys.phases.at(phase_name).comps.at(comp_name).k_well;
          }
      }
  }

  void
  calculate_locals()
  {
    calculate_omega_phase();
    calculate_sum_sq_eta();
    calculate_h();
    calculate_dhdeta();
    calculate_local_diffusivity();
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
                                                   comp.dmudt.grad * userInputs.dtValue);
        var_index++;
      }
    for (auto &[phase_name, op] : op_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                op.eta.val +
                                                  op.detadt.val * userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   op.detadt.grad * userInputs.dtValue);
        var_index++;
      }
  }

  // void
  // submit_pp_fields(
  //   variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
  //   uint                                                            &var_index)
  //{}
};

#endif