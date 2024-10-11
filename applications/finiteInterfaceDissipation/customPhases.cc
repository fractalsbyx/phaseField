// This file will be created or modified by AMMBER
#include "SystemContainer.h"

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_A : public PhaseFieldContainer<dim, degree> {
public:
    const FieldContainer<dim>& x_CU;
    const FieldContainer<dim>& x_SI;

    FieldContainer<dim>& dfdx_CU;
    FieldContainer<dim>& dfdx_SI;

public:
    Phase_A(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list),
          x_CU(this->comp_data["CU"].x_data),
          x_SI(this->comp_data["SI"].x_data),
          dfdx_CU(this->comp_data["CU"].dfdx),
          dfdx_SI(this->comp_data["SI"].dfdx)
    {}

    void calculate_G() override {
        this->phase_free_energy = (x_CU.val * log(x_CU.val) + x_SI.val * log(x_SI.val)) + 2.5 * x_CU.val * x_SI.val;
    }

    void calculate_dfdx() override {
        dfdx_CU.val = (log(x_CU.val) + 1.0) + 2.5 * x_SI.val;
        dfdx_SI.val = (log(x_SI.val) + 1.0) + 2.5 * x_CU.val;

        dfdx_CU.grad = (1.0/(x_CU.val))*x_CU.grad + (2.5)*x_SI.grad;
        dfdx_SI.grad = (2.5)*x_CU.grad + (1.0/(x_SI.val))*x_SI.grad;
    }
};

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class Phase_B : public PhaseFieldContainer<dim, degree> {
public:
    const FieldContainer<dim>& x_CU;
    const FieldContainer<dim>& x_SI;

    FieldContainer<dim>& dfdx_CU;
    FieldContainer<dim>& dfdx_SI;

public:
    Phase_B(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list),
          x_CU(this->comp_data["CU"].x_data),
          x_SI(this->comp_data["SI"].x_data),
          dfdx_CU(this->comp_data["CU"].dfdx),
          dfdx_SI(this->comp_data["SI"].dfdx)
    {}

    void calculate_G() override {
        this->phase_free_energy = (x_CU.val * log(x_CU.val) + x_SI.val * log(x_SI.val)) + 1.5 * x_CU.val * x_SI.val +6.0;
    }

    void calculate_dfdx() override {
        dfdx_CU.val = (log(x_CU.val) + 1.0) + 1.5 * x_SI.val;
        dfdx_SI.val = (log(x_SI.val) + 1.0) + 1.5 * x_CU.val;

        dfdx_CU.grad = (1.0/(x_CU.val))*x_CU.grad + (1.5)*x_SI.grad;
        dfdx_SI.grad = (1.5)*x_CU.grad + (1.0/(x_SI.val))*x_SI.grad;
    }
};

// 
template <int dim, int degree>
SystemContainer<dim, degree>::SystemContainer(const IsothermalSystem& _isoSys, variableContainer<dim,degree,scalarValue>& _variable_list) :
        isoSys(_isoSys), variable_list(_variable_list){
    // For all phase names
    phase_fields.insert({"Phase_A", new Phase_A<dim,degree>(isoSys, "Phase_A", phase_fields, variable_list)});
    phase_fields.insert({"Phase_B", new Phase_B<dim,degree>(isoSys, "Phase_B", phase_fields, variable_list)});
    //phase_fields.insert({"FCC", new FCC<dim,degree>(isoSys, "FCC", phase_fields, variable_list)});
    //phase_fields.insert({"BCC", new BCC<dim,degree>(isoSys, "BCC", phase_fields, variable_list)});
    //phase_fields.insert({"HCP", new HCP<dim,degree>(isoSys, "HCP", phase_fields, variable_list)});
    // etc.
}
