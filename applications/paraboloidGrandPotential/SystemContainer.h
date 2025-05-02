#ifndef SYSTEMCONTAINER_H
#define SYSTEMCONTAINER_H

#include <deal.II/base/exceptions.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "FieldContainer.h"
#include "ParaboloidSystem.h"

#include <core/userInputParameters.h>
#include <core/variableContainer.h>
#include <map>
#include <string>

// template <typename number>
// using boost_vector = boost::numeric::ublas::vector<number>;
// template <typename number>
// using boost_symmet = boost::numeric::ublas::symmetric_matrix<number>;
#include "LinAlg.h"

namespace boost::numeric::ublas
{

  template <int dim>
  struct promote_traits<dealii::VectorizedArray<double>,
                        dealii::Tensor<1, dim, dealii::VectorizedArray<double>>>
  {
    using promote_type = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
  };

  template <int dim>
  struct promote_traits<dealii::Tensor<1, dim, dealii::VectorizedArray<double>>,
                        dealii::VectorizedArray<double>>
  {
    using promote_type = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
  };

} // namespace boost::numeric::ublas

template <int dim, int degree>
class SystemContainer
{
public:
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;

  /**
   * @brief Data structure to hold the phase data
   */
  struct PhaseData
  {
    FieldContainer<dim>               omega;
    FieldContainer<dim>               h;
    boost_vector<FieldContainer<dim>> delta_mu;
    boost_vector<FieldContainer<dim>> delta_c;
    boost_vector<FieldContainer<dim>> c_phase;
  };

  /**
   * @brief Data structure to hold the order parameter data
   */
  struct OPData
  {
    FieldContainer<dim>              eta;
    Variation<dim>                   detadt;
    std::vector<FieldContainer<dim>> dhdeta;
  };

  /**
   * @brief Pointer to the system parameters
   */
  const ParaboloidSystem &isoSys;
  /**
   * @brief Pointer to the PRISMS-PF parameters
   */
  const userInputParameters<dim> &userInputs;

  /**
   * @brief Values associated with each phase
   */
  std::vector<PhaseData> phase_data;

  boost_vector<FieldContainer<dim>> mu;
  boost_vector<Variation<dim>>      dmudt;
  boost_symmet<scalarValue>         M;

  /**
   * @brief Values associated with each order parameter
   */
  std::vector<std::pair<uint, OPData>> op_data;
  /**
   * @brief Sum of squares of the order parameters
   */
  FieldContainer<dim> sum_sq_eta;

  /**
   * Constructor
   */
  SystemContainer(const ParaboloidSystem &sys, const userInputParameters<dim> &inputs)
    : isoSys(sys)
    , userInputs(inputs)
    , phase_data(std::vector<PhaseData>(isoSys.phases.size()))
    , mu(boost_vector<FieldContainer<dim>>(isoSys.num_comps))
    , dmudt(boost_vector<Variation<dim>>(isoSys.num_comps))
    , M(boost_symmet<scalarValue>(isoSys.num_comps))
    , op_data({})
    , sum_sq_eta({})
  {}

  /**
   * @brief Initialize the fields for the PDE
   * @param variable_list The variable list
   * @param var_index The starting index for the block of fields
   */
  void
  initialize_fields_explicit(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys.order_params.size());
    mu.resize(isoSys.num_comps);
    for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
      {
        mu[comp_index].val  = variable_list.get_scalar_value(var_index);
        mu[comp_index].grad = variable_list.get_scalar_gradient(var_index);
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

  /**
   * @brief Initialize the fields needed for the postprocess
   * @param variable_list The variable list
   * @param var_index The starting index for the block of fields
   */
  void
  initialize_fields_postprocess(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys.order_params.size());
    mu.resize(isoSys.num_comps);
    for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
      {
        mu[comp_index].val = variable_list.get_scalar_value(var_index);
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
  calculate_deltas()
  {
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases[phase_index];
        PhaseData                     &phase      = phase_data[phase_index];
        phase.delta_mu                            = mu - phase_info.B_well;
        phase.delta_c = prod(phase_info.suscept, phase.delta_mu);
        phase.c_phase = phase_info.c_ref + phase.delta_c;
      }
  }

  /**
   * @brief Calculate the grand potential density for each phase
   */
  void
  calculate_omega_phase()
  {
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases[phase_index];
        PhaseData                     &phase      = phase_data[phase_index];
        phase.omega =
          -0.5 * inner_prod(phase.delta_mu, prod(phase_info.suscept, phase.delta_mu)) -
          inner_prod(phase_info.c_ref, mu) + phase_info.D_well;
      }
  }

  /**
   * @brief Calculate the sum of squares of the order parameters
   */
  void
  calculate_sum_sq_eta()
  {
    sum_sq_eta.val = dealii::make_vectorized_array(0.);
    for (const auto &[phase_index, op] : op_data)
      {
        sum_sq_eta += op.eta * op.eta;
      }
  }

  /**
   * @brief Calculate the phase fraction, h, for each phase
   * @details h_i = eta_i^2 / sum(eta^2)
   */
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

  /**
   * @brief Calculate the derivative of the phase fraction, h, with respect to eta
   */
  void
  calculate_dhdeta()
  {
    for (auto &[alpha_index, op] : op_data)
      {
        for (uint beta_index = 0; beta_index < op.dhdeta.size(); beta_index++)
          {
            PhaseData           &beta   = phase_data[beta_index];
            FieldContainer<dim> &dhdeta = op.dhdeta[beta_index];
            dhdeta.val                  = dealii::make_vectorized_array(0.);
            if (alpha_index == beta_index)
              {
                dhdeta += 2.0 * op.eta;
              }
            dhdeta -= 2.0 * op.eta * beta.h;
            dhdeta /= sum_sq_eta;
          }
      }
  }

  /**
   * @brief Calculate the time evolution of the order parameters
   */
  void
  calculate_detadt()
  {
    for (auto &[alpha_index, op] : op_data)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases.at(alpha_index);

        double m     = 6.00 * phase_info.sigma / isoSys.l_int;
        double kappa = 0.75 * phase_info.sigma * isoSys.l_int;
        double L     = 4.00 * phase_info.mu_int / isoSys.l_int / 3.00;

        // Interface term
        Variation<dim> interface_term;
        interface_term.val =
          m * (op.eta.val * op.eta.val * op.eta.val - op.eta.val +
               2. * 1.5 * op.eta.val * (sum_sq_eta.val - op.eta.val * op.eta.val));
        interface_term.vec = -kappa * op.eta.grad;

        // Chemical term
        // This is a variation, but has no vector term.
        scalarValue chemical_term = dealii::make_vectorized_array(0.);
        for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
          {
            const PhaseData &beta = phase_data.at(beta_index);
            chemical_term += beta.omega.val * op.dhdeta.at(beta_index).val;
          }
        op.detadt = -L * (interface_term + chemical_term);
      }
  }

  /**
   * @brief Calculate the local chemical mobility
   */
  void
  calculate_local_mobility()
  {
    M = boost_symmet<scalarValue>(isoSys.num_comps);
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        PhaseData                     &phase      = phase_data[phase_index];
        const ParaboloidSystem::Phase &phase_info = isoSys.phases.at(phase_index);
        const double                  &D          = phase_info.D;
        const boost_symmet<double>    &suscept    = phase_info.suscept;
        /* M += phase.h.val * isoSys.phases.at(phase_index).D *
             isoSys.phases.at(phase_index).suscept; */ // NOTATION
        for (uint comp_row = 0; comp_row < isoSys.num_comps; comp_row++)
          {
            for (uint comp_col = 0; comp_col <= comp_row; comp_col++)
              {
                M(comp_row, comp_col) += phase.h.val * D * suscept(comp_row, comp_col);
              }
          }
      }
  }

  /**
   * @brief Calculate the time evolution of the chemical potential
   * @details Calculates the evolution of the composition, then converts to chemical
   * potential through the susceptibility
   */
  void
  calculate_dmudt()
  {
    // Get the local susceptibility
    boost_symmet<FieldContainer<dim>> local_suscept(isoSys.num_comps, isoSys.num_comps);
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys.phases[phase_index];
        // local_suscept_inv += phase_info.suscept * phase_data[phase_index].h; //
        // NOTATION
        for (uint comp_row = 0; comp_row < isoSys.num_comps; comp_row++)
          {
            for (uint comp_col = 0; comp_col <= comp_row; comp_col++)
              {
                local_suscept(comp_row, comp_col) +=
                  phase_info.suscept(comp_row, comp_col) * phase_data[phase_index].h;
              }
          }
      }
    boost_matrix<FieldContainer<dim>> local_suscept_inv =
      inverse_from_cholesky(local_suscept);

    // // Get just the gradient of mu
    // boost_vector<scalarGrad> mu_grad(isoSys.num_comps);
    // for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
    //   {
    //     mu_grad[comp_index] = mu[comp_index].grad;
    //   }

    // Flux term
    // dmudt += prod(M, mu_grad); // NOTATION
    // boost_vector<scalarGrad> flux = prod(M, mu_grad);
    // for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
    //  {
    //    dmudt[comp_index] += flux[comp_index];
    //  }
    dmudt = boost_vector<Variation<dim>>(isoSys.num_comps);
    for (uint comp_row = 0; comp_row < isoSys.num_comps; comp_row++)
      {
        for (uint comp_col = 0; comp_col < isoSys.num_comps; comp_col++)
          {
            dmudt[comp_row].vec += M(comp_row, comp_col) * mu[comp_col].grad;
          }
      }

    // Partition/conservation term
    for (auto &[phase_index, op] : op_data)
      {
        boost_vector<FieldContainer<dim>> dcdeta_sum(isoSys.num_comps);
        for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
          {
            dcdeta_sum += op.dhdeta.at(beta_index) * (phase_data[beta_index].c_phase);
          }
        // dmudt -= op.detadt * dcdeta_sum; // NOTATION
        for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
          {
            dmudt[comp_index] -= op.detadt * dcdeta_sum[comp_index];
          }
      }
    // Convert from dcdt to dmudt
    dmudt = prod(local_suscept_inv, dmudt);
  }

  /**
   * @brief Calculate the information needed to solve the evolution equations in the
   * proper order
   */
  void
  calculate_locals()
  {
    calculate_deltas();
    calculate_omega_phase();
    calculate_sum_sq_eta();
    calculate_h();
    calculate_dhdeta();
    calculate_local_mobility();
  }

  /**
   * @brief Submit the fields to PRISMS-PF
   * @param variable_list The PRISMS-PF variable list
   * @param var_index The starting index for the block of fields
   */
  void
  submit_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                            &var_index)
  {
    for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                mu[comp_index].val +
                                                  dmudt[comp_index].val *
                                                    userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -dmudt[comp_index].vec *
                                                     userInputs.dtValue);
        var_index++;
      }
    for (auto &[phase_index, op] : op_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                op.eta.val +
                                                  op.detadt.val * userInputs.dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -op.detadt.vec * userInputs.dtValue);
        var_index++;
      }
  }

  /**
   * @brief Submit the post-processed fields to PRISMS-PF
   * @param pp_variable_list The PRISMS-PF variable list
   * @param pp_index The starting index for the block of fields
   */
  void
  submit_pp_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
    uint                                                            &pp_index)
  {
    boost_vector<FieldContainer<dim>> c_loc(isoSys.num_comps);
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const PhaseData &phase = phase_data.at(phase_index);
        c_loc += phase.h * phase.c_phase;
      }
    for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
      {
        pp_variable_list.set_scalar_value_term_RHS(pp_index, c_loc[comp_index].val);
        pp_index++;
      }
  }
};

#endif