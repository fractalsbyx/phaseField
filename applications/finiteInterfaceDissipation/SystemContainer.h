#include "PhaseFieldContainer.h"
#include <map>


template <int dim, int degree>
class SystemContainer{
public:
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    const IsothermalSystem& isoSys;
    variableContainer<dim,degree,scalarValue>& variable_list;
    std::map<std::string, PhaseFieldContainer<dim, degree>*> phase_fields;
    SystemContainer(const IsothermalSystem& _isoSys,
                    variableContainer<dim,degree,scalarValue>& _variable_list) :
                    isoSys(_isoSys), variable_list(_variable_list) {
    }
    ~SystemContainer(){
        for(auto& [key, phase_field] : phase_fields){
            delete phase_field;
        }
    }

    void initialize_fields(){
        uint var_index = 0;
        for(auto& [key, phase_field] : phase_fields){
            phase_field.initialize_fields(var_index);
        }
    }

    void solve(){
        for(auto& [key, phase_field] : phase_fields){
            phase_field.solve(phase_fields);
        }
    }

    void submit_fields(const double& dt){
        uint var_index = 0;
        for(auto& [key, phase_field] : phase_fields){
            phase_field.submit_fields(var_index, dt);
        }
    }
};