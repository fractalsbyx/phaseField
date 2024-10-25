// This file will be created or modified by AMMBER
#include "SystemContainer.h"

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_A : public PhaseFieldContainer<dim, degree> {
public:
    Phase_A(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list)
    {}

    inline void calculate_free_energy() override {
        const FieldContainer<dim>& x_CU = this->comp_data["CU"].x_data;
        const FieldContainer<dim>& x_SI = this->comp_data["SI"].x_data;

        FieldContainer<dim>& dfdx_CU = this->comp_data["CU"].dfdx;
        FieldContainer<dim>& dfdx_SI = this->comp_data["SI"].dfdx;

        this->phase_free_energy = (x_CU.val * x_CU.val) + 0.5 * x_SI.val;

        dfdx_CU.val = (2.0*x_CU.val);
        dfdx_SI.val = 0.5;

        dfdx_CU.grad = (2.0)*x_CU.grad + (0.0)*x_SI.grad;
        dfdx_SI.grad = (0.0)*x_CU.grad + (0.0)*x_SI.grad;
    }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_B : public PhaseFieldContainer<dim, degree> {
public:
    Phase_B(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list)
    {}

    inline void calculate_free_energy() override {
        const FieldContainer<dim>& x_CU = this->comp_data["CU"].x_data;
        const FieldContainer<dim>& x_SI = this->comp_data["SI"].x_data;

        FieldContainer<dim>& dfdx_CU = this->comp_data["CU"].dfdx;
        FieldContainer<dim>& dfdx_SI = this->comp_data["SI"].dfdx;

        this->phase_free_energy = (x_SI.val * x_SI.val) + 0.5 * x_CU.val;

        dfdx_CU.val = 0.5;
        dfdx_SI.val = (2.0*x_SI.val);

        dfdx_CU.grad = (0.0)*x_CU.grad + (0.0)*x_SI.grad;
        dfdx_SI.grad = (0.0)*x_CU.grad + (2.0)*x_SI.grad;
    }
};

// 
template <int dim, int degree>
inline SystemContainer<dim, degree>::SystemContainer(const IsothermalSystem& _isoSys, variableContainer<dim,degree,scalarValue>& _variable_list) :
        isoSys(_isoSys), variable_list(_variable_list){
    // For all phase names
    phase_fields.insert({"Phase_A", new Phase_A<dim,degree>(isoSys, "Phase_A", phase_fields, variable_list)});
    phase_fields.insert({"Phase_B", new Phase_B<dim,degree>(isoSys, "Phase_B", phase_fields, variable_list)});
}
template class SystemContainer<2,1>;
template class SystemContainer<2,2>;
template class SystemContainer<2,3>;
template class SystemContainer<3,1>;
template class SystemContainer<3,2>;
template class SystemContainer<3,3>;