#ifndef PHASEFIELDCONTAINER_H
#define PHASEFIELDCONTAINER_H

#include <deal.II/base/vectorization.h>

#include "FieldContainer.h"
#include "IsothermalSystem.h"

#include "../../include/variableContainer.h"
#include <map> //
#include <set>
#include <string> //

template <int dim, int degree>
class SystemContainer;

template <int dim>
struct CompData
{
  FieldContainer<dim> x_data;
  FieldContainer<dim> dfdx;
  FieldContainer<dim> dxdt;
};

constexpr double PI = 3.141592653589793238;

template <int dim, int degree>
class PhaseFieldContainer
{
public:
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)
  double epsilon = 1.0e-10;

  PhaseFieldContainer(
    const IsothermalSystem                                          &isoSys,
    const std::string                                               &phase_name,
    const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : isoSys(isoSys)
    , info(isoSys.phases.at(phase_name))
    , phase_fields(phase_fields)
  {}

  virtual ~PhaseFieldContainer() = default;

  // defined in customPhases.cc
  virtual void
  calculate_chem_energy()
  {}

  void
  initialize_fields(uint                                              &var_index,
                    const variableContainer<dim, degree, scalarValue> &variable_list)
  {
    // Phase Value
    phi.val  = variable_list.get_scalar_value(var_index);
    phi.grad = variable_list.get_scalar_gradient(var_index);
    var_index++;
    // Components Value
    for (const auto &[comp_name, comp_info] : info.comps)
      {
        comp_data[comp_name].x_data.val  = variable_list.get_scalar_value(var_index);
        comp_data[comp_name].x_data.grad = variable_list.get_scalar_gradient(var_index);
        var_index++;
      }
  }

  void
  calculate_dfdphi(const scalarValue &sum_sq_phi)
  {
    dfdphi.val = info.m0 * ((phi.val * phi.val * phi.val - phi.val) +
                            3. * phi.val * (sum_sq_phi - phi.val * phi.val)) +
                 2. * phi.val * phase_free_energy;
    dfdphi.grad = info.kappa * -phi.grad;
  }

  void
  calculate_dphidt(const FieldContainer<dim> &backflux_term)
  {
    dphidt = info.mu * (-dfdphi - field_x_variation(phi, backflux_term));
  }

  void
  calculate_dxdt()
  {
    FieldContainer<dim> sum_M_x, sum_minus_M_dfdx;
    sum_minus_M_dfdx *= 0.;
    sum_M_x *= 0.;
    for (auto &[i, i_alpha] : comp_data)
      {
        sum_M_x += info.comps.at(i).M * i_alpha.x_data;
        sum_minus_M_dfdx -= info.comps.at(i).M * i_alpha.dfdx;
      }
    FieldContainer<dim> backflux_term =
      field_x_field(sum_M_x.inverse(), sum_minus_M_dfdx);
    for (auto &[i, i_alpha] : comp_data)
      {
        i_alpha.dxdt *= 0.0;
        // Spatial flux
        i_alpha.dxdt.grad +=
          info.comps.at(i).M *
          (-i_alpha.dfdx.grad - (i_alpha.x_data.val * backflux_term.grad));
        i_alpha.dxdt.val -=
          info.comps.at(i).M *
          (-i_alpha.dfdx.grad - (i_alpha.x_data.val * backflux_term.grad)) * -phi.grad /
          (phi.val + epsilon);

        // Internal relaxation (eq. 16)
        scalarValue         pairsum1 = constV(0.0);
        FieldContainer<dim> pairsum2;
        pairsum2.val = constV(0.0);
        pairsum2.grad *= 0.0;
        for (const auto &[beta_name, beta] : phase_fields)
          {
            if (beta != this) // Technically unnecessary
              {
                const CompData<dim> &i_beta = beta->comp_data.at(i);
                pairsum1 += beta->phi.val * (i_beta.dfdx.val - i_alpha.dfdx.val);
                pairsum2.val +=
                  (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.val;
                pairsum2.grad +=
                  (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.grad;
                pairsum2.val -=
                  (i_beta.x_data.grad - i_alpha.x_data.grad) * beta->dphidt.grad;
              }
          }
        i_alpha.dxdt.val += isoSys.comp_info.at(i).P * pairsum1 - pairsum2.val;
        i_alpha.dxdt.grad += -pairsum2.grad;
      }
  }

  void
  submit_fields(uint                                        &var_index,
                variableContainer<dim, degree, scalarValue> &variable_list,
                const double                                &dt)
  {
    variable_list.set_scalar_value_term_RHS(var_index, phi.val + dt * dphidt.val);
    variable_list.set_scalar_gradient_term_RHS(var_index, -dt * dphidt.grad);
    var_index++;
    for (auto &[i, i_data] : comp_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                i_data.x_data.val + dt * i_data.dxdt.val);
        variable_list.set_scalar_gradient_term_RHS(var_index, -dt * i_data.dxdt.grad);
        var_index++;
      }
  }

  // protected:
  //  References to phase object
  const IsothermalSystem                                          &isoSys;
  const Phase                                                     &info;
  const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields;

  // Object containing values for components of this phase
  std::map<std::string, CompData<dim>> comp_data;
  std::set<std::string>                comp_names;
  // Free energy G for this phase at its composition
  scalarValue phase_free_energy;

  FieldContainer<dim> phi;
  FieldContainer<dim> dfdphi;
  FieldContainer<dim> dphidt;
};

#endif