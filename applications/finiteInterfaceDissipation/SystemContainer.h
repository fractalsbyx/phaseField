#include "PhaseFieldContainer.h"
#include "IsothermalSystem.h"
#include <set>


template <int dim, int degree>
class SystemContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    SystemContainer(const IsothermalSystem& isoSys, variableContainer<dim,degree,scalarValue>& variable_list){
        uint var_index = 0;
        for(const Phase& phase : isoSys.phases){
            phase_containers.insert(PhaseFieldContainer<dim, degree>(phase, variable_list, var_index));
        }
    }
    std::set<PhaseFieldContainer<dim, degree>> phase_containers;
};