#include "Phase.h"
#include <deal.II/base/vectorization.h>
//#include "../../include/variableContainer.h"


template <int dim>
struct FieldContainer{
    dealii::VectorizedArray<double> val;
    dealii::Tensor<1, dim, dealii::VectorizedArray<double>> grad;
};

template <int dim, int degree>
class PhaseFieldContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    PhaseFieldContainer(const Phase& phase, variableContainer<dim,degree,scalarValue>& variable_list, uint& var_index) :
                                    name(phase.name),
                                    comps(phase.comps),
                                    variable_list(variable_list),
                                    first_var_index(var_index){
        // Phase Value
        
        // Components Value
        for (const auto& comp : comps){
            x_data[comp].val = variable_list.get_scalar_value(var_index);
            x_data[comp].grad = variable_list.get_scalar_gradient(var_index++);
        }
        phi.val *= 0.0;
        phi.grad *= 0.0;
    }
    ~PhaseFieldContainer(){}

    public:

    // Custom START
    void calculate_dfdx(){
        for(const auto& comp : comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }
    void calculate_G(){
        for(const auto& comp : comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }
    // Custom END

    void calculate_dfdphi1(){
        for(const auto& phi : order_parameters){
            dfdphi[comp].val = constV(0.0);
            dfdphi[comp].grad *= constV(0.0);
        }
    }
    void calculate_dfdphi2(){
        for(auto& op : order_parameters){
            op.dfdphi.grad *= constV(0.0);
            op.dfdphi.val = op.dfdphi.val 
        }
    }

    void calculate_dxdt(){
        scalarGrad temp;
        temp *= 0.0;
        const scalarValue VAL_ZERO = constV(0.0);
        const scalarGrad VEC_ZERO = temp;

        for(const auto& i : comps){
            // Spatial flux
            for(const auto& i : comps){
                dxdt[i].grad += Vm*Vm * M[i][j] * dfdx[i].grad;
                dxdt[i].val += Vm*Vm * M[i][j] * phi.grad * dfdx.grad / phi.val;
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

    void submit_dxdt(const double& dt){
        uint var_index = first_var_index;
        for (const std::string& comp : phase.comps){
            variable_list.set_scalar_value_term_RHS(var_index, y_data[s][constituent].val + dt * dydt_data[s][constituent].val);
            variable_list.set_scalar_value_term_RHS(var_index++,                          - dt * dydt_data[s][constituent].grad);
        }
    }


    protected:
    // References to phase object
    std::string& name;
    std::unordered_set<std::string>& comps;

    variableContainer<dim,degree,scalarValue>& variable_list;
    const uint first_var_index;

    std::unordered_map<std::string, FieldContainer<dim>> x_data;
    std::unordered_map<std::string, FieldContainer<dim>> dfdx_data;
    std::unordered_map<std::string, FieldContainer<dim>> dxdt_data;

    FieldContainer<dim> phi;
    FieldContainer<dim> dphidt;
};
