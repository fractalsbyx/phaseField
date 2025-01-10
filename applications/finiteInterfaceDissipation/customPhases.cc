// This file will be created or modified by AMMBER
#include "SystemContainer.h"

constexpr double k_b = 8.617333e-5; // eV/K
constexpr double T   = 1775.;       // K
constexpr double Va  = 0.01;        // nm3

constexpr double kT = k_b * T / Va;

constexpr double omega = 90.; // eV/nm3

constexpr double L_CU = 11.5;  // eV/nm3
constexpr double L_TI = 11.8;  // eV/nm3
constexpr double L_AA = 17.6;  // eV/nm3
constexpr double T_CU = 1358.; // K
constexpr double T_TI = 1941.; // K
constexpr double T_AA = 3290.; // K
constexpr double K_CU = L_CU * (T - T_CU) / T_CU;
constexpr double K_TI = L_TI * (T - T_TI) / T_TI;
constexpr double K_AA = L_AA * (T - T_AA) / T_AA;

template <int dim, int degree>
dealii::VectorizedArray<double>
PhaseFieldContainer<dim, degree>::M_ij(const std::string &i, const std::string &j)
{
  return info.comps.at(i).M * comp_data.at(i).x_data.val *
         (double(i == j) - comp_data.at(j).x_data.val); // 0.5 * (info.comps.at(i).M
                                                        // + info.comps.at(j).M);
                                                        // // fast bad approximation
}

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class SOLID : public PhaseFieldContainer<dim, degree>
{
public:
  SOLID(const IsothermalSystem                                          &isoSys,
        const std::string                                               &phase_name,
        const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_CU    = this->comp_data["CU"].x_data;
    FieldContainer<dim>       &dfdx_CU = this->comp_data["CU"].dfdx;

    const FieldContainer<dim> &x_TI    = this->comp_data["TI"].x_data;
    FieldContainer<dim>       &dfdx_TI = this->comp_data["TI"].dfdx;

    this->phase_free_energy =
      kT * (                                                                     //
             x_CU.val * std::log(x_CU.val) +                                     //
             x_TI.val * std::log(x_TI.val) +                                     //
             ((1.0 - x_CU.val - x_TI.val) * std::log(1.0 - x_CU.val - x_TI.val)) //
             ) +                                                                 //
      omega * x_CU.val * (1.0 - x_CU.val - x_TI.val) +                           //
      K_CU * x_CU.val +
      K_TI * x_TI.val + K_AA * (1.0 - x_CU.val - x_TI.val);

    dfdx_CU.val = kT * (std::log(x_CU.val) -                         //
                        std::log(1.0 - x_CU.val - x_TI.val)) +       //
                  omega * ((1.0 - x_CU.val - x_TI.val) - x_CU.val) + //
                  K_CU -
                  K_AA;
    dfdx_TI.val = kT * (std::log(x_TI.val) -                   //
                        std::log(1.0 - x_CU.val - x_TI.val)) + //
                  omega * (-x_CU.val) +                        //
                  K_TI -
                  K_AA;

    dfdx_CU.grad =
      (kT / x_CU.val + kT / (1.0 - x_CU.val - x_TI.val) - omega * 2.0) * x_CU.grad + //
      (kT / (1.0 - x_CU.val - x_TI.val) - omega) * x_TI.grad;
    dfdx_TI.grad = (kT / (1.0 - x_CU.val - x_TI.val) - omega) * x_CU.grad + //
                   (kT / x_TI.val + kT / (1.0 - x_CU.val - x_TI.val)) * x_TI.grad;

    // this->volumetrize_free_energy();
    this->nondimensionalize_free_energy();
  }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class LIQUID : public PhaseFieldContainer<dim, degree>
{
public:
  LIQUID(const IsothermalSystem                                          &isoSys,
         const std::string                                               &phase_name,
         const std::map<std::string, PhaseFieldContainer<dim, degree> *> &phase_fields)
    : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields)
  {}

  inline void
  calculate_free_energy() override
  {
    const FieldContainer<dim> &x_CU    = this->comp_data["CU"].x_data;
    FieldContainer<dim>       &dfdx_CU = this->comp_data["CU"].dfdx;

    const FieldContainer<dim> &x_TI    = this->comp_data["TI"].x_data;
    FieldContainer<dim>       &dfdx_TI = this->comp_data["TI"].dfdx;

    this->phase_free_energy =
      kT * (                                                                     //
             x_CU.val * std::log(x_CU.val) +                                     //
             x_TI.val * std::log(x_TI.val) +                                     //
             ((1.0 - x_CU.val - x_TI.val) * std::log(1.0 - x_CU.val - x_TI.val)) //
             ) +                                                                 //
      omega * x_CU.val * (1.0 - x_CU.val - x_TI.val);

    dfdx_CU.val = kT * (std::log(x_CU.val) -                   //
                        std::log(1.0 - x_CU.val - x_TI.val)) + //
                  omega * ((1.0 - x_CU.val - x_TI.val) - x_CU.val);
    dfdx_TI.val = kT * (std::log(x_TI.val) -                   //
                        std::log(1.0 - x_CU.val - x_TI.val)) + //
                  omega * (-x_CU.val);

    dfdx_CU.grad =
      (kT / x_CU.val + kT / (1.0 - x_CU.val - x_TI.val) - omega * 2.0) * x_CU.grad + //
      (kT / (1.0 - x_CU.val - x_TI.val) - omega) * x_TI.grad;
    dfdx_TI.grad = (kT / (1.0 - x_CU.val - x_TI.val) - omega) * x_CU.grad + //
                   (kT / x_TI.val + kT / (1.0 - x_CU.val - x_TI.val)) * x_TI.grad;

    // this->volumetrize_free_energy();
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
  phase_fields.insert({"SOLID", new SOLID<dim, degree>(isoSys, "SOLID", phase_fields)});
  phase_fields.insert(
    {"LIQUID", new LIQUID<dim, degree>(isoSys, "LIQUID", phase_fields)});
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