// This file will be created or modified by AMMBER
#include "SystemContainer.h"

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class ALPHA : public PhaseFieldContainer<dim, degree> {
public:
    ALPHA(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list)
    {}
    #define constV(a) dealii::make_vectorized_array(a)
    inline void calculate_free_energy() override {
        const FieldContainer<dim>& x_AA = this->comp_data["AA"].x_data;
        const FieldContainer<dim>& x_BB = this->comp_data["BB"].x_data;

        FieldContainer<dim>& dfdx_AA = this->comp_data["AA"].dfdx;
        FieldContainer<dim>& dfdx_BB = this->comp_data["BB"].dfdx;

        this->phase_free_energy = (x_BB.val * x_BB.val);

        dfdx_AA.val = constV(0.);
        dfdx_BB.val = 2.0*x_BB.val;

        dfdx_AA.grad = (0.0)*x_AA.grad + (0.0)*x_BB.grad;
        dfdx_BB.grad = (0.0)*x_AA.grad + (2.0)*x_BB.grad;
    }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class LIQUID : public PhaseFieldContainer<dim, degree> {
public:
    LIQUID(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list)
    {}
    #define constV(a) dealii::make_vectorized_array(a)
    inline void calculate_free_energy() override {
        const FieldContainer<dim>& x_AA = this->comp_data["AA"].x_data;
        const FieldContainer<dim>& x_BB = this->comp_data["BB"].x_data;

        FieldContainer<dim>& dfdx_AA = this->comp_data["AA"].dfdx;
        FieldContainer<dim>& dfdx_BB = this->comp_data["BB"].dfdx;

        this->phase_free_energy = (x_AA.val * x_AA.val) + constV(0.1);

        dfdx_AA.val = 2.0*x_AA.val;
        dfdx_BB.val = constV(0.);

        dfdx_AA.grad = (2.0)*x_AA.grad + (0.0)*x_BB.grad;
        dfdx_BB.grad = (0.0)*x_AA.grad + (0.0)*x_BB.grad;
    }
};

// 
template <int dim, int degree>
inline SystemContainer<dim, degree>::SystemContainer(const IsothermalSystem& _isoSys,
                                                    variableContainer<dim,degree,scalarValue>& _variable_list,
                                                    const userInputParameters<dim>& _userInputs) :
                                                    isoSys(_isoSys),
                                                    variable_list(_variable_list),
                                                    userInputs(_userInputs){
    // For all phase names
    phase_fields.insert({"ALPHA", new ALPHA<dim,degree>(isoSys, "ALPHA", phase_fields, variable_list)});
    phase_fields.insert({"LIQUID", new LIQUID<dim,degree>(isoSys, "LIQUID", phase_fields, variable_list)});
}
template class SystemContainer<2,1>;
template class SystemContainer<2,2>;
template class SystemContainer<2,3>;
template class SystemContainer<3,1>;
template class SystemContainer<3,2>;
template class SystemContainer<3,3>;