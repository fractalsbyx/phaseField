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
double epsilon = 1.0e-10;

template <int dim, int degree>
class PhaseFieldContainer{
public:
    using scalarValue = dealii::VectorizedArray<double>;
    using scalarGrad =  dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
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
    // defined in customPhases.cc
    virtual void calculate_dfdx(){}
    virtual void calculate_G(){}

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
    // Equation 33
    void calculate_I(){
        I.val = phi.val*PI*PI/isoSys.eta;
        I.grad = phi.grad;
    }

    void calculate_dxdt(){
        for(auto& [i, i_alpha] : comp_data){
            i_alpha.dxdt.val = constV(0.0);
            i_alpha.dxdt.grad *= 0.0;
            // Spatial flux
            for(auto& [j, j_alpha] : comp_data){
                i_alpha.dxdt.grad += isoSys.Vm*isoSys.Vm * M_ij(i,j) * j_alpha.dfdx.grad;
                i_alpha.dxdt.val -= -isoSys.Vm*isoSys.Vm * M_ij(i,j) * j_alpha.dfdx.grad * phi.grad / abs(phi.val+epsilon);
            }
            // Internal relaxation (eq. 16)
            scalarValue pairsum1 = constV(0.0);
            FieldContainer<dim> pairsum2;
            pairsum2.val = constV(0.0);
            pairsum2.grad *= 0.0;
            for (const auto& [beta_name, beta] : phase_fields){
                if (beta != this){ // Technically unnecessary 
                    const CompData<dim>& i_beta = beta->comp_data.at(i);
                    pairsum1 += beta->phi.val * (i_beta.dfdx.val - i_alpha.dfdx.val);
                    pairsum2.val += (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.val;
                    pairsum2.grad += (i_beta.x_data.val - i_alpha.x_data.val) * beta->dphidt.grad;
                    pairsum2.val -= (i_beta.x_data.grad - i_alpha.x_data.grad) * beta->dphidt.grad;
                }
            }
            i_alpha.dxdt.val += isoSys.comp_info.at(i).P * pairsum1 - pairsum2.val;
            i_alpha.dxdt.grad += - pairsum2.grad;
        }
    }
    // Equation 37
    FieldContainer<dim> K(const PhaseFieldContainer& beta){
        double mu_ab = mu(*this, beta); //TODO
        FieldContainer<dim> symmetric_term;
        symmetric_term.val = 4.0*(double)isoSys.N*isoSys.eta*(phi.val+beta.phi.val);
        symmetric_term.grad = 4.0*(double)isoSys.N*isoSys.eta*(phi.grad+beta.phi.grad);
        FieldContainer<dim> denom_sum_term;
        denom_sum_term.val = constV(0.0);
        denom_sum_term.grad *= 0.0;
        for (auto& [i, i_alpha] : comp_data){
            const CompData<dim>& i_beta = beta.comp_data.at(i);
            denom_sum_term.val +=   (i_alpha.x_data.val - i_beta.x_data.val)*
                                    (i_alpha.x_data.val - i_beta.x_data.val)/
                                    isoSys.comp_info.at(i).P;
            denom_sum_term.grad +=  2.0 * 
                                    (i_alpha.x_data.val - i_beta.x_data.val)*
                                    (i_alpha.x_data.grad - i_beta.x_data.grad)/
                                    isoSys.comp_info.at(i).P;
        }
        FieldContainer<dim> denominator;
        denominator.val = symmetric_term.val + (mu_ab * PI*PI * denom_sum_term.val) /*+ epsilon */;
        denominator.grad = symmetric_term.grad + (mu_ab * PI*PI * denom_sum_term.grad);
        FieldContainer<dim> K_out;
        K_out.val = mu_ab*symmetric_term.val/denominator.val;
        K_out.grad = mu_ab*((symmetric_term.grad * denominator.val) - 
                                (symmetric_term.val * denominator.grad))/
                                (denominator.val * denominator.val);
        //std::cout << "Sym: " << symmetric_term << ", N: " << isoSys.N << ", musum: " << (mu_ab * PI*PI * denom_sum_term) << ", ";
        return K_out;
    }
    // Equation 38
    scalarValue delta_G_phi_ab(const PhaseFieldContainer& beta){
        scalarValue sum_term = constV(0.0);
        for (const auto& [i, i_alpha] : comp_data){
            const CompData<dim>& i_beta = beta.comp_data.at(i);
            sum_term += (phi.val*i_alpha.dfdx.val + beta.phi.val*i_beta.dfdx.val)*
                        (i_beta.x_data.val - i_alpha.x_data.val);
        }
        sum_term /= phi.val + beta.phi.val /*+ epsilon*/;
        //std::cout << "dG: " << beta.phase_free_energy - phase_free_energy - sum_term << ", ";
        return beta.phase_free_energy - phase_free_energy - sum_term;
    }
    // Equation 39
    void calculate_dphidt(){
        dphidt.val = constV(0.0);
        dphidt.grad *= 0.0;
        //std::cout << "ZV: " << dphidt.val << ", ";
        //std::cout << "ZG: " << dphidt.grad << ", ";
        for(const auto& [beta_name, beta] : phase_fields){
            if (beta != this){
                FieldContainer<dim> inner_sum_term;
                inner_sum_term.val = constV(0.0);
                inner_sum_term.grad *= 0.0;
                for(const auto& [gamma_name, gamma] : phase_fields){
                    if(gamma != this && gamma != beta){
                        inner_sum_term.val += (sigma(*beta, *gamma) - sigma(*this, *gamma)) * gamma->I.val;
                        inner_sum_term.grad += (sigma(*beta, *gamma) - sigma(*this, *gamma)) * gamma->I.grad;
                    }
                }
                FieldContainer<dim> K_ab = K(*beta);
                dphidt.val += K_ab.val * 
                            (sigma(*this, *beta)*(I.val - beta->I.val) +
                            inner_sum_term.val + 0.25*PI*PI*delta_G_phi_ab(*beta)/isoSys.eta);
                //std::cout << "dphiV: " << dphidt.val << ", ";
                dphidt.grad += K_ab.val * 
                            (sigma(*this, *beta)*(I.grad - beta->I.grad) +
                            inner_sum_term.grad);
                dphidt.val -= K_ab.grad * 
                            (sigma(*this, *beta)*(I.grad - beta->I.grad) +
                            inner_sum_term.grad);
                //std::cout << "dphiG: " << dphidt.grad << ", ";
            }
        }
        dphidt.val /= (double)isoSys.N;
        dphidt.grad /= (double)isoSys.N;
    }

    inline double sigma(const PhaseFieldContainer<dim, degree>& alpha, const PhaseFieldContainer<dim, degree>& beta){
        return 0.5*(alpha.info.sigma + beta.info.sigma);
    }
    inline double M_ij(const std::string& i, const std::string& j){
        return 0.5*(info.comps.at(i).M + info.comps.at(j).M); // fast bad approximation
    }
    inline double mu(const PhaseFieldContainer<dim, degree>& alpha, const PhaseFieldContainer<dim, degree>& beta){
        return 0.5*(alpha.info.mu + beta.info.mu); // fast bad approximation
    }

    void calculate_locals(){
        calculate_G();
        calculate_dfdx();
        calculate_I();
    }

    void solve(){
        calculate_dphidt();
        calculate_dxdt();
    }

    void submit_fields(uint& var_index, const double& dt){
        //std::cout << "dt: " << dt << ", ";
        //std::cout << "Va: " << dphidt.val << ", ";
        //std::cout << "Gr: " << dphidt.grad << ", ";
        //std::cout << "DphiV: " << phi.val + dt * dphidt.val << ", ";
        //std::cout << "DphiG: " << - dt * dphidt.grad << ", ";
        variable_list.set_scalar_value_term_RHS(var_index, phi.val + dt * dphidt.val);
        variable_list.set_scalar_gradient_term_RHS(var_index++,    - dt * dphidt.grad);
        for (auto& [i, i_data] : comp_data){
            variable_list.set_scalar_value_term_RHS(var_index, i_data.x_data.val + dt * i_data.dxdt.val);
            variable_list.set_scalar_gradient_term_RHS(var_index++,              - dt * i_data.dxdt.grad);
        }
    }

protected:
    // References to phase object
    const IsothermalSystem& isoSys;
    const Phase& info;
    const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields;
    variableContainer<dim,degree,scalarValue>& variable_list;

    std::map<std::string, CompData<dim>> comp_data;
    std::set<std::string> comp_names;

    scalarValue phase_free_energy;

    FieldContainer<dim> phi;
    FieldContainer<dim> dphidt;
    FieldContainer<dim> I;
};
