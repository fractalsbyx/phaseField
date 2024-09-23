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
        x_data = std::unordered_map<std::string, FieldContainer<dim>>;
        dfdx_data = std::unordered_map<std::string, FieldContainer<dim>>;
        dxdt_data = std::unordered_map<std::string, FieldContainer<dim>>;
        for (const auto& comp : comps){
            x_data[comp].val = variable_list.get_scalar_value(var_index);
            x_data[comp].grad = variable_list.get_scalar_gradient(var_index++);
        }
        phi_phase.val *= 0.0;
        phi_phase.grad *= 0.0;
    }
    ~PhaseFieldContainer(){}

    public:

    void calculate_dfdx(){
        for(const auto& comp : comps){
            dfdx[comp].val = constV(0.0);
            dfdx[comp].grad *= constV(0.0);
        }
    }

    void calculate_dxdt(){
        scalarGrad temp;
        temp *= 0.0;
        const scalarValue VAL_ZERO = constV(0.0);
        const scalarGrad VEC_ZERO = temp;

        for(const auto& comp : comps){
            dxdt[comp].grad += M[comp]*dfdx[comp].grad;
            dxdt[comp].val += phi_phase.grad * dfdx.grad * [comp]/phi_phase.val;

            dxdt[comp].val += P[comp] *;
        }
        return;
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

    FieldContainer<dim> phi_phase;
};
