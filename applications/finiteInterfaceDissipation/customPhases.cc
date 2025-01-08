// This file will be created or modified by AMMBER
#include "SystemContainer.h"

constexpr double omega = 0.0;

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_A : public PhaseFieldContainer<dim, degree>
{
public:
  Phase_A(const IsothermalSystem                                          &isoSys,
          const std::string                                               &phase_name,
          const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_CU = this->comp_data["CU"].x_data;

    FieldContainer<dim> &dfdx_CU = this->comp_data["CU"].dfdx;

    /* this->phase_free_energy = (x_CU.val * x_CU.val) + 0.5 * (constV(1.0) - x_CU.val);

    dfdx_CU.val = 2.0 * (x_CU.val - constV(0.25));

    dfdx_CU.grad = (2.0) * x_CU.grad; */

    this->phase_free_energy = x_CU.val * std::log(x_CU.val)                   //
                              + ((1.0 - x_CU.val) * std::log(1.0 - x_CU.val)) //
                              + omega * x_CU.val * (1.0 - x_CU.val)           //
                              + (x_CU.val);

    dfdx_CU.val = std::log(x_CU.val) - std::log(1.0 - x_CU.val) //
                  + omega * ((1.0 - x_CU.val) - x_CU.val)       //
                  + 1.0;

    dfdx_CU.grad = (1.0 / x_CU.val + 1.0 / (1.0 - x_CU.val) + omega * -2.0) * x_CU.grad;

    this->volumetrize_free_energy();
    this->nondimensionalize_free_energy();
  }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_B : public PhaseFieldContainer<dim, degree>
{
public:
  Phase_B(const IsothermalSystem                                          &isoSys,
          const std::string                                               &phase_name,
          const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_CU = this->comp_data["CU"].x_data;

    FieldContainer<dim> &dfdx_CU = this->comp_data["CU"].dfdx;

    /* this->phase_free_energy =
      ((constV(1.0) - x_CU.val) * (constV(1.0) - x_CU.val)) + 0.5 * x_CU.val;

    dfdx_CU.val = 2.0 * (x_CU.val - constV(0.75));

    dfdx_CU.grad = (2.0) * x_CU.grad; */

    this->phase_free_energy = x_CU.val * std::log(x_CU.val)                   //
                              + ((1.0 - x_CU.val) * std::log(1.0 - x_CU.val)) //
                              + omega * x_CU.val * (1.0 - x_CU.val)           //
                              + (1.0 - x_CU.val);

    dfdx_CU.val = std::log(x_CU.val) - std::log(1.0 - x_CU.val) //
                  + omega * ((1.0 - x_CU.val) - x_CU.val)       //
                  - 1.0;

    dfdx_CU.grad = (1.0 / x_CU.val + 1.0 / (1.0 - x_CU.val) + omega * -2.0) * x_CU.grad;

    this->volumetrize_free_energy();
    this->nondimensionalize_free_energy();
  }
};

//
template <int dim, int degree>
inline SystemContainer<dim, degree>::SystemContainer(
  const IsothermalSystem         &_isoSys,
  const userInputParameters<dim> &_userInputs)
  : isoSys(_isoSys)
  , userInputs(_userInputs)
{
  // For all phase names
  phase_fields.insert(
    {"Phase_A", new Phase_A<dim, degree>(isoSys, "Phase_A", phase_fields)});
  phase_fields.insert(
    {"Phase_B", new Phase_B<dim, degree>(isoSys, "Phase_B", phase_fields)});
}
template class SystemContainer<2, 1>;
template class SystemContainer<2, 2>;
template class SystemContainer<2, 3>;
template class SystemContainer<2, 4>;
template class SystemContainer<2, 5>;
template class SystemContainer<2, 6>;
template class SystemContainer<3, 1>;
template class SystemContainer<3, 2>;
template class SystemContainer<3, 3>;
template class SystemContainer<3, 4>;
template class SystemContainer<3, 5>;
template class SystemContainer<3, 6>;