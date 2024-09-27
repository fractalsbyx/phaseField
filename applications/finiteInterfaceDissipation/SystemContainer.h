#include "PhaseFieldContainer.h"
#include "IsothermalSystem.h"
#include <map>


template <int dim, int degree>
class SystemContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    SystemContainer(const IsothermalSystem& isoSys, variableContainer<dim,degree,scalarValue>& variable_list){};
    ~SystemContainer(){
        for(PhaseFieldContainer<dim, degree>* [key, phase_field] : phase_fields){
            delete phase_field;
        }
    }
    std::map<std::string, PhaseFieldContainer<dim, degree>*> phase_fields;
};