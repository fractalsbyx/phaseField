#include "IsothermalSystem.h"
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
public:
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    PhaseFieldContainer(const IsothermalSystem& isoSys, const std::string& phase_name,
                        const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                        variableContainer<dim,degree,scalarValue>& variable_list) :
                                    isoSys(isoSys),
                                    info(isoSys.phases.at(phase_name)),
                                    phase_fields(phase_fields),
                                    variable_list(variable_list){
    }
    virtual ~PhaseFieldContainer(){}

    void initialize_fields(uint& var_index){
        // Phase Value
        phi.val = variable_list.get_scalar_value(var_index);
        phi.grad = variable_list.get_scalar_gradient(var_index++);
        // Components Value
        for (const auto& [comp_name, comp_info] : info.comps){
            comp_data[comp_name].x_data.val = variable_list.get_scalar_value(var_index);
            comp_data[comp_name].x_data.grad = variable_list.get_scalar_gradient(var_index++);
        }
    }

    // Custom START
    virtual void calculate_dfdx(){
    }
    virtual void calculate_G(){
    }
    // Custom END

    void calculate_dxdt(){
        scalarGrad temp;
        temp *= 0.0;
        const scalarValue VAL_ZERO = constV(0.0);
        const scalarGrad VEC_ZERO = temp;

        for(auto& [i, i_alpha] : comp_data){
            // Spatial flux
            for(auto& [j, j_alpha] : comp_data){
                i_alpha.dxdt.grad += isoSys.Vm*isoSys.Vm * M_ij(i,j) * j_alpha.dfdx.grad;
                i_alpha.dxdt.val +=  isoSys.Vm*isoSys.Vm * M_ij(i,j) * phi.grad * j_alpha.dfdx.grad / phi.val;
            }
            // Internal relaxation (eq. 16)
            scalarValue pairsum1 = constV(0.0);
            FieldContainer<dim> pairsum2;
            pairsum2.val = constV(0.0);
            pairsum2.grad *= 0.0;
            for (const auto& [beta_name, beta] : phase_fields){
                const CompData<dim>& i_beta = beta->comp_data.at(i);
                pairsum1 += beta->phi.val * (i_beta.dfdx.val - i_alpha.dfdx.val);
                pairsum2.val += phi.val * (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.val;
                pairsum2.val += (phi.val * (i_beta.x_data.grad - i_alpha.x_data.grad)
                                +phi.grad * (i_beta.x_data.val - i_alpha.x_data.val)
                                ) * beta->dphidt.grad;
                pairsum2.grad -= phi.val * (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.grad;
            }
            i_alpha.dxdt.val += isoSys.comp_info.at(i).P * phi.val * pairsum1 + pairsum2.val;
            i_alpha.dxdt.grad += pairsum2.grad;
        }
    }
    // Equation 37
    scalarValue K_ab(const PhaseFieldContainer& beta){
        scalarValue mu_ab; //TODO
        scalarValue symmetric_term = 4.0*(double)isoSys.N*isoSys.eta*(phi.val+beta.phi.val);
        scalarValue denom_sum_term = constV(0.0);
        for (auto& [i, i_alpha] : comp_data){
            const CompData<dim>& i_beta = beta.comp_data.at(i);
            denom_sum_term +=   (i_alpha.x_data.val - i_beta.x_data.val)*
                                (i_alpha.x_data.val - i_beta.x_data.val)/
                                isoSys.comp_info.at(i).P;
        }
        scalarValue denominator = 1.0 + (mu_ab * PI*PI * denom_sum_term)/symmetric_term;
        return mu_ab/denominator;
    }
    // Equation 38
    scalarValue delta_G_phi_ab(const PhaseFieldContainer& beta){
        scalarValue sum_term = constV(0.0);
        for (auto& [i, i_alpha] : comp_data){
            const CompData<dim>& i_beta = beta.comp_data.at(i);
            sum_term += (phi.val*i_alpha.dfdx.val + beta.phi.val*i_beta.dfdx.val)*
                        (i_beta.x_data.val - i_alpha.x_data.val);
        }
        sum_term /= phi.val + beta.phi.val;
        return beta.phase_free_energy - phase_free_energy - sum_term;
    }
    // Equation 39
    void calculate_dphidt(){
        dphidt.val = constV(0.0);
        dphidt.grad *= 0.0;
        for(const auto& [beta_name, beta] : phase_fields){
            FieldContainer<dim> inner_sum_term;
            inner_sum_term.val = constV(0.0);
            inner_sum_term.grad *= 0.0;
            for(const auto& [gamma_name, gamma] : phase_fields){
                if(gamma != this && gamma != beta){
                    inner_sum_term.val += (sigma(*beta, *gamma) - sigma(*this, *gamma)) * gamma->I.val;
                    inner_sum_term.grad += (sigma(*beta, *gamma) - sigma(*this, *gamma)) * gamma->I.grad;
                }
            }
            dphidt.val += K_ab(*beta) * 
                        (sigma(*this, *beta)*(I.val - beta->I.val) +
                        inner_sum_term.val);
            dphidt.grad += K_ab(*beta) * 
                        (sigma(*this, *beta)*(I.grad - beta->I.grad) +
                        inner_sum_term.grad);
                        
        }
        dphidt.val /= (double)isoSys.N;
        dphidt.grad /= (double)isoSys.N;
    }

    void calculate_I(){
        I.val = phi.val*PI*PI/isoSys.eta;
        I.grad = phi.grad;
    }

    inline double sigma(const PhaseFieldContainer<dim, degree>& alpha, const PhaseFieldContainer<dim, degree>& beta){
        return 0.5*(alpha.info.sigma + beta.info.sigma);
    }
    inline double M_ij(const std::string& i, const std::string& j){
        return 0.5*(info.comps.at(i).M + info.comps.at(j).M); // fast bad approximation
    }

    void solve(){
        calculate_G();
        calculate_dfdx();
        calculate_I();
        calculate_dphidt();
        calculate_dxdt();
    }

    void submit_fields(uint& var_index, const double& dt){
        variable_list.set_scalar_value_term_RHS(var_index, phi.val + dt * dphidt.val);
        variable_list.set_scalar_gradient_term_RHS(var_index++,    - dt * dphidt.grad);
        for (auto& [i, i_data] : comp_data){
            variable_list.set_scalar_value_term_RHS(var_index, i_data.x_data.val + dt * i_data.dxdt.val);
            variable_list.set_scalar_gradient_term_RHS(var_index++,              - dt * i_data.dxdt.grad);
        }
    }

protected:
    // References to phase object
    const Phase& info;
    const IsothermalSystem& isoSys;
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
    FieldContainer<dim> I;
};
