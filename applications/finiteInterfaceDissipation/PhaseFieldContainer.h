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
    FieldContainer<dim> dfdx;
    FieldContainer<dim> dxdt;
};
constexpr double PI = 3.141592653589793238;

template <int dim, int degree>
class PhaseFieldContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    PhaseFieldContainer(const Phase& phase, std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields, variableContainer<dim,degree,scalarValue>& variable_list) :
                                    info(phase),
                                    phase_fields(phase_fields),
                                    variable_list(variable_list){
    }
    virtual ~PhaseFieldContainer(){}

public:

    void initialize_fields(uint& var_index){
        // Phase Value
        phi.val = variable_list.get_scalar_value(var_index);
        phi.grad = variable_list.get_scalar_gradient(var_index++);
        // Components Value
        for (const auto& comp : info.comps){
            comp_data[comp].x_data.val = variable_list.get_scalar_value(var_index);
            comp_data[comp].x_data.grad = variable_list.get_scalar_gradient(var_index++);
        }
    }

    // Custom START
    virtual void calculate_dfdx(){
        for(const auto& comp : info.comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }
    virtual void calculate_G(){
        for(const auto& comp : info.comps){
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

        for(auto& [i, i_alpha] : comp_data){
            // Spatial flux
            for(const auto& i : info.comps){
                i_alpha.dxdt.grad += Vm*Vm * M_ij * i_alpha.dfdx.grad;
                i_alpha.dxdt.val +=  Vm*Vm * M_ij * phi.grad * i_alpha.dfdx.grad / phi.val;
            }
            // Internal relaxation (eq. 16)
            scalarValue pairsum1 = constV(0.0);
            FieldContainer pairsum2 = constV(0.0);
            for (const auto& beta : phaselist){
                const CompData& i_beta = beta.comp_data.at(i);
                pairsum1 += beta.phi.val * (i_beta.dfdx.val - i_alpha.dfdx.val);
                pairsum2.val += phi.val * (i_beta.x_data.val - i_alpha.x_data.val) * beta.dphidt.val;
                pairsum2.val += (phi.val * (i_beta.x_data.grad - i_alpha.x_data.grad)
                                +phi.grad * (i_beta.x_data.val - i_alpha.x_data.val)
                                ) * beta.dphidt.grad;
                pairsum2.grad -= phi.val * (i_beta.x_data.val - i_alpha.x_data.val) * beta.dphidt.grad;
            }
            i_alpha.dxdt.val += P[i] * phi.val * pairsum1 + pairsum2.val;
            i_alpha.dxdt.grad += pairsum2.grad;
        }
    }
    // Equation 37
    scalarValue K_ab(const PhaseFieldContainer& beta){
        scalarValue mu_ab; //TODO
        scalarValue symmetric_term = 4.0*isoSys.N*isoSys.eta*(phi.val+beta.phi.val);
        scalarValue denom_sum_term = constV(0.0);
        for (auto& [i, i_alpha] : comp_data){
            const CompData& i_beta = beta.comp_data.at(i);
            denom_sum_term +=   (i_alpha.x_data.val - i_beta.x_data.val)*
                                (i_alpha.x_data.val - i_beta.x_data.val)/
                                P[i];
        }
        scalarValue denominator = 1.0 + (mu_ab * PI*PI * denom_sum_term)/symmetric_term;
        return mu_ab/denominator;
    }
    // Equation 38
    scalarValue delta_G_phi_ab(const PhaseFieldContainer& beta){
        scalarValue sum_term = constV(0.0);
        for (auto& [i, i_alpha] : comp_data){
            const CompData& i_beta = beta.comp_data.at(i);
            sum_term += (phi.val*i_alpha.dfdx.val + beta.phi.val*i_beta.dfdx.val)*
                        (i_beta.x_data.val - i_alpha.x_data.val);
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
                    inner_sum_term += (sigma(beta, gamma) - sigma(*this, gamma)) * gamma.I.val;
                }
            }
            dphidt.val += K_ab(beta) * 
                        (sigma(*this, beta)*(I.val - beta.I.val) +
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
            variable_list.set_scalar_value_term_RHS(var_index, x_data[comp].val + dt * dxdt[comp].val);
            variable_list.set_scalar_value_term_RHS(var_index++,                - dt * dxdt[comp].grad);
        }
    }

protected:
    // References to phase object
    const Phase& info;
    const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields;
    variableContainer<dim,degree,scalarValue>& variable_list;

    std::unordered_map<std::string, CompData<dim>> comp_data;
    std::set<std::string> comp_names;
    // std::unordered_map<std::string, FieldContainer<dim>> x_data;
    // std::unordered_map<std::string, FieldContainer<dim>> dfdx;
    // std::unordered_map<std::string, FieldContainer<dim>> dxdt;

    scalarValue phase_free_energy;

    FieldContainer<dim> phi;
    FieldContainer<dim> dphidt;
};
