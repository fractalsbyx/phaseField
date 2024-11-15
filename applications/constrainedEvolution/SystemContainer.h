#ifndef SYSTEMCONTAINER_H
#define SYSTEMCONTAINER_H

#include "FieldContainer.h"
#include "PhaseFieldContainer.h"

#include "../../include/userInputParameters.h"
#include <map>    //
#include <string> //

template <int dim, int degree>
class SystemContainer
{
public:
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)

  // Members
  const IsothermalSystem                                   &isoSys;
  const userInputParameters<dim>                           &userInputs;
  std::map<std::string, PhaseFieldContainer<dim, degree> *> phase_fields;

  FieldContainer<dim> sum_sq_phi, sum_mu_phi, sum_minus_mu_dfdphi;

  SystemContainer(const IsothermalSystem         &_isoSys,
                  const userInputParameters<dim> &_userInputs);

  ~SystemContainer()
  {
    for (auto &[key, phase_field] : phase_fields)
      {
        delete phase_field;
      }
  }

  void
  initialize_fields(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list)
  {
    // Zero the uninitialized globals
    sum_sq_phi          = 0. * sum_sq_phi;
    sum_mu_phi          = sum_sq_phi;
    sum_minus_mu_dfdphi = sum_sq_phi;

    // Get solution fields
    uint var_index = 0;
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->initialize_fields(var_index, variable_list);
        sum_sq_phi += field_x_field(phase_field->phi, phase_field->phi);
      }
  }

  void
  calculate_locals()
  {
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_chem_energy();
      }
    sum_mu_phi.val = constV(0.);
    sum_mu_phi.grad *= 0.;
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_dfdphi(sum_sq_phi.val);
      }
  }

  void
  solve()
  {
    for (auto &[key, phase_field] : phase_fields)
      {
        sum_mu_phi += phase_field->info.mu * phase_field->phi;
        sum_minus_mu_dfdphi -= phase_field->info.mu * phase_field->dfdphi;
      }
    FieldContainer<dim> backflux_term =
      field_x_variation(sum_mu_phi.inverse(), sum_minus_mu_dfdphi);
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_dphidt(backflux_term);
      }
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_dxdt();
      }
  }

  void
  submit_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list)
  {
    uint var_index = 0;
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->submit_fields(var_index, variable_list, userInputs.dtValue);
      }
  }

  void
  submit_pp_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list)
  {
    uint        var_index = 0;
    scalarValue pp_comp;
    for (const auto &[comp_name, other] : isoSys.comp_info)
      {
        pp_comp = constV(0.);
        for (auto &[key, phase_field] : phase_fields)
          {
            pp_comp +=
              phase_field->phi.val * phase_field->comp_data[comp_name].x_data.val;
          }
        pp_variable_list.set_scalar_value_term_RHS(var_index, pp_comp);
        ++var_index;
      }
  }
};

#endif