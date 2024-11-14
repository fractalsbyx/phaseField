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
  FieldContainer<dim>                                       sum_sq_psi;

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
    sum_sq_psi.val = constV(0.);
    sum_sq_psi.grad *= 0.;
    uint var_index = 0;
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->initialize_fields(var_index, variable_list);
        sum_sq_psi += phase_field->phi;
      }
    scalarValue magnitude = sqrt(sum_sq_psi.val);
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->psi.val /= magnitude;
        phase_field->psi.grad /= magnitude;
        phase_field->phi.val /= sum_sq_psi.val;
        phase_field->phi.grad /= sum_sq_psi.val;
      }
    magnitude /= magnitude;
    sum_sq_psi.val = constV(1.);
    sum_sq_psi.grad *= 0.;
  }

  void
  calculate_locals()
  {
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_locals();
      }
  }

  void
  solve()
  {
    FieldContainer<dim> constraint_term;
    constraint_term.val *= 0.;
    constraint_term.grad *= 0.;
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_dfdpsi(sum_sq_psi.val);
        constraint_term += field_x_variation(phase_field->psi,
                                             -phase_field->info.mu * phase_field->dfdpsi);
      }
    for (auto &[key, phase_field] : phase_fields)
      {
        phase_field->calculate_dpsidt(constraint_term);
        // phase_field->calculate_dpsidt(FieldContainer<dim> {constV(0.), 0. *
        // constraint_term.grad});
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
              phase_field->psi.val * phase_field->comp_data[comp_name].x_data.val;
          }
        pp_variable_list.set_scalar_value_term_RHS(var_index, pp_comp);
        ++var_index;
      }
  }
};

#endif