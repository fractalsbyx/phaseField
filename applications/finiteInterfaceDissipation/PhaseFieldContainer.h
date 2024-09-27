#include "Phase.h"
#include <deal.II/base/vectorization.h>
#include "../../include/variableContainer.h"


template <int dim>
struct FieldContainer{
    dealii::VectorizedArray<double> val;
    dealii::Tensor<1, dim, dealii::VectorizedArray<double>> grad;
};

template <int dim>
struct CompData{
    FieldContainer<dim> x_data;
    FieldContainer<dim> dfdx_data;
    FieldContainer<dim> dxdt_data;
};
constexpr double PI = 3.141592653589793238;

template <int dim, int degree>
class PhaseFieldContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    PhaseFieldContainer(const Phase& phase, variableContainer<dim,degree,scalarValue>& variable_list, uint& var_index) :
                                    info(phase),
                                    variable_list(variable_list),
                                    first_var_index(var_index){
        // Phase Value
        phi.val = variable_list.get_scalar_value(var_index);
        phi.grad = variable_list.get_scalar_gradient(var_index++);
        // Components Value
        for (const auto& comp : comps){
            x_data[comp].val = variable_list.get_scalar_value(var_index);
            x_data[comp].grad = variable_list.get_scalar_gradient(var_index++);
        }
    }
    virtual ~PhaseFieldContainer(){}

public:

    // Custom START
    virtual void calculate_dfdx(){
        for(const auto& comp : comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }
    virtual void calculate_G(){
        for(const auto& comp : comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }
    // Custom END

    void calculate_dxdt(){
        scalarGrad temp;
        temp *= 0.0;
        const scalarValue VAL_ZERO = constV(0.0);
        const scalarGrad VEC_ZERO = temp;

        for(const auto& i : comps){
            // Spatial flux
            for(const auto& i : comps){
                dxdt[i].grad += Vm*Vm * M_ij * dfdx[i].grad;
                dxdt[i].val += Vm*Vm * M_ij * phi.grad * dfdx.grad / phi.val;
            }
            // Internal relaxation (eq. 16)
            scalarValue pairsum1 = constV(0.0);
            FieldContainer pairsum2 = constV(0.0);
            for (const auto& beta : phaselist){
                pairsum1 += beta.phi.val * (beta.dfdx_data[i].val - dfdx_data[i].val);
                pairsum2.val += phi.val * (beta.x_data[i].val - x_data[i].val) * beta.dphidt.val;
                pairsum2.val += (phi.val * (beta.x_data[i].grad - x_data[i].grad)
                                +phi.grad * (beta.x_data[i].val - x_data[i].val)
                                ) * beta.dphidt.grad;
                pairsum2.grad -= phi.val * (beta.x_data[i].val - x_data[i].val) * beta.dphidt.grad;
            }
            dxdt[i].val += P[i] * phi.val * pairsum1 + pairsum2.val;
            dxdt[i].grad += pairsum2.grad;
        }
    }
    // Equation 37
    scalarValue K_ab(const PhaseFieldContainer& beta){
        scalarValue mu_ab;
        scalarValue symmetric_term = 4.0*isoSys.N*isoSys.eta*(phi.val+beta.phi.val);
        scalarValue denom_sum_term = constV(0.0);
        for (const auto& comp : comps){
            denom_sum_term +=   (x_data[comp].val - beta.x_data[comp].val)*
                                (x_data[comp].val - beta.x_data[comp].val)/
                                P[comp];
        }
        scalarValue denominator = 1.0 + (mu_ab * PI*PI * denom_sum_term)/symmetric_term;
        return mu_ab/denominator;
    }
    // Equation 38
    scalarValue delta_G_phi_ab(const PhaseFieldContainer& beta){
        scalarValue sum_term = constV(0.0);
        for (const auto& comp : comps){
            sum_term += (phi.val*dfdx[comp].val + beta.phi.val*beta.dfdx[comp].val)*
                        (beta.x_data[comp].val - x_data[comp].val);
        }
        sum_term /= phi.val + beta.phi.val;
        return beta.phase_free_energy - phase_free_energy - sum_term;
    }
    // Equation 39
    void calculate_dphidt(){
        dphidt.val = constV(0.0);
        for(const PhaseFieldContainer& beta : phaselist){
            scalarValue inner_sum_term = constV(0.0);
            for(const PhaseFieldContainer& gamma : phaselist){
                if(&gamma != this && &gamma != &beta){
                    inner_sum_term += (sigma_beta_gamma - sigma_alpha_gamma) * gamma.I.val;
                }
            }
            dphidt.val += K_ab(beta) * 
                        (sigma_ab*(I.val - beta.I.val) +
                        inner_sum_term);
        }
        dphidt.val /= isoSys.N;
        dphidt.grad /= isoSys.N;
    }

    inline double sigma(PhaseFieldContainer<dim, degree>& alpha, PhaseFieldContainer<dim, degree>& beta){
        return 0.5*(alpha.info.sigma + beta.info.sigma);
    }

    void solve(){
        calculate_G();
        calculate_dfdx();
        calculate_dphidt();
        calculate_dxdt();
    }

    void submit_fields(uint& var_index){
        variable_list.set_scalar_value_term_RHS(var_index, phi.val + dt * dphidt.val);
        variable_list.set_scalar_value_term_RHS(var_index++,       - dt * dphidt.grad);
        for (const std::string& comp : phase.comps){
            variable_list.set_scalar_value_term_RHS(var_index, x_data[comp].val + dt * dxdt_data[comp].val);
            variable_list.set_scalar_value_term_RHS(var_index++,                - dt * dxdt_data[comp].grad);
        }
    }

    protected:
    // References to phase object
    const Phase& info;

    variableContainer<dim,degree,scalarValue>& variable_list;
    const uint first_var_index;

    std::unordered_map<std::string, CompData<dim>> comps;
    std::set<std::string> comp_names;
    std::unordered_map<std::string, FieldContainer<dim>> x_data;
    std::unordered_map<std::string, FieldContainer<dim>> dfdx_data;
    std::unordered_map<std::string, FieldContainer<dim>> dxdt_data;

    scalarValue phase_free_energy;

    FieldContainer<dim> phi;
    FieldContainer<dim> dphidt;
};
