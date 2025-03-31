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
    FieldContainer<dim>              eta;
    FieldContainer<dim>              detadt;
    std::vector<FieldContainer<dim>> dhdeta;
  };

  const ParaboloidSystem         &isoSys;
  const userInputParameters<dim> &userInputs;

  std::vector<PhaseData>               phase_data;
  std::vector<CompData>                comp_data;
  std::vector<std::pair<uint, OPData>> op_data;
  FieldContainer<dim>                  sum_sq_eta;

  /**
   * Constructor
   */
  SystemContainer(const ParaboloidSystem &sys, const userInputParameters<dim> &inputs)
    : isoSys(sys)
    , userInputs(inputs)
    , phase_data(std::vector<PhaseData>(isoSys.phases.size()))
    , comp_data(std::vector<CompData>(isoSys.comp_names.size()))
    , op_data({})
    , sum_sq_eta({})
  {}

  ~SystemContainer()
  {}

  void
  initialize_fields_explicit(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys.order_params.size());
    for (uint comp_index = 0; comp_index < isoSys.comp_names.size(); comp_index++)
      {
        comp_data[comp_index].mu.val  = variable_list.get_scalar_value(var_index);
        comp_data[comp_index].mu.grad = variable_list.get_scalar_gradient(var_index);
        var_index++;
      }
    for (const auto &phase_index : isoSys.order_params)
      {
        OPData op;
        op.eta.val  = variable_list.get_scalar_value(var_index);
        op.eta.grad = variable_list.get_scalar_gradient(var_index);
        op.dhdeta.resize(isoSys.phases.size());
        op_data.push_back({phase_index, op});
        var_index++;
      }
  }

  void
  initialize_fields_postprocess(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys.order_params.size());
    for (uint comp_index = 0; comp_index < isoSys.comp_names.size(); comp_index++)
      {
        comp_data[comp_index].mu.val = variable_list.get_scalar_value(var_index);
        var_index++;
      }
    for (const auto &phase_index : isoSys.order_params)
      {
        OPData op;
        op.eta.val = variable_list.get_scalar_value(var_index);
        op.dhdeta.resize(isoSys.phases.size());
        op_data.push_back({phase_index, op});
        var_index++;
      }
  }

  void
  calculate_omega_phase()
  {
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases[phase_index];
        PhaseData                     &phase      = phase_data[phase_index];
        phase.omega.val                           = phase_info.f_min;
        for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
          {
            const CompData                        &comp = comp_data.at(comp_index);
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              phase_info.comps.at(comp_index);
            phase.omega +=
              -comp.mu * comp.mu / (2.0 * isoSys.Vm * isoSys.Vm * comp_info.k_well) -
              comp.mu * comp_info.c_min / isoSys.Vm;
          }
      }
  }

  void
  calculate_sum_sq_eta()
  {
    sum_sq_eta.val = constV(0.);
    for (const auto &[phase_index, op] : op_data)
      {
        sum_sq_eta += op.eta * op.eta;
      }
  }

  void
  calculate_h()
  {
    for (auto &[phase_index, op] : op_data)
      {
        phase_data[phase_index].h += op.eta * op.eta;
      }
    for (auto &phase : phase_data)
      {
        phase.h /= sum_sq_eta;
      }
  }

  void
  calculate_dhdeta()
  {
    for (auto &[alpha_index, op] : op_data)
      {
        for (uint beta_index = 0; beta_index < op.dhdeta.size(); beta_index++)
          {
            PhaseData           &beta   = phase_data[beta_index];
            FieldContainer<dim> &dhdeta = op.dhdeta[beta_index];
            dhdeta.val                  = constV(0.);
            if (alpha_index == beta_index)
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
    for (auto &[alpha_index, op] : op_data)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases.at(alpha_index);

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
        for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
          {
            const PhaseData &beta = phase_data.at(beta_index);
            // NOTE: Only multiply values
            chemical_term += beta.omega.val * op.dhdeta.at(beta_index).val;
          }
        op.detadt = -L * (interface_term + chemical_term);
      }
  }

  void
  calculate_local_mobility()
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData &comp = comp_data[comp_index];
        comp.M         = constV(0.);
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            PhaseData &phase = phase_data[phase_index];
            comp.M += isoSys.phases.at(phase_index).D * phase.h.val /
                      (isoSys.Vm * isoSys.Vm *
                       isoSys.phases.at(phase_index).comps.at(comp_index).k_well);
          }
      }
  }

  void
  calculate_dmudt()
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData           &comp = comp_data[comp_index];
        FieldContainer<dim> chi_AA;
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            PhaseData                             &phase = phase_data[phase_index];
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys.phases.at(phase_index).comps.at(comp_index);
            chi_AA += phase.h / comp_info.k_well / isoSys.Vm / isoSys.Vm;
          }
        comp.dmudt.val = constV(0.);
        // Flux term
        comp.dmudt.grad = -comp.M * -comp.mu.grad;

        // Partitioning term
        for (auto &[phase_index, op] : op_data)
          {
            FieldContainer<dim> drhodeta_sum;
            for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
              {
                auto &comp_info = isoSys.phases.at(beta_index).comps.at(comp_index);
                drhodeta_sum +=
                  op.dhdeta.at(beta_index) *
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
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData &comp = comp_data[comp_index];
        variable_list.set_scalar_value_term_RHS(var_index,
                                                comp.mu.val +
                                                  comp.dmudt.val * userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -comp.dmudt.grad * userInputs.dtValue);
        var_index++;
      }
    for (auto &[phase_index, op] : op_data)
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
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData   &comp = comp_data[comp_index];
        scalarValue c    = constV(0.);
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            const PhaseData                       &phase = phase_data.at(phase_index);
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys.phases.at(phase_index).comps.at(comp_index);
            c += phase.h.val * comp_info.c_min;
            c += phase.h.val * comp.mu.val / comp_info.k_well / isoSys.Vm;
          }
        pp_variable_list.set_scalar_value_term_RHS(pp_index, c);
        pp_index++;
      }
  }
};

#endif