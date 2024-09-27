// This file will be created or modified by AMMBER
#include "PhaseFieldContainer.h"

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class ExamplePhase : public PhaseFieldContainer<dim, degree>{
private:
    const scalarValue& x_CU = x_data["CU"].val;
    const scalarValue& x_SI = x_data["SI"].val;
    const scalarValue& x_MG = x_data["MG"].val;

    scalarValue& dfdx_CU = dfdx_data["CU"].val;
    scalarValue& dfdx_SI = dfdx_data["SI"].val;
    scalarValue& dfdx_MG = dfdx_data["MG"].val;
public:
    phase_name = "ExamplePhase";
    ExamplePhase() :    x_CU(x_data["CU"].val),
                        x_SI(x_data["SI"].val),
                        x_MG(x_data["MG"].val),
                        dfdx_CU(dfdx_data["CU"].val),
                        dfdx_SI(dfdx_data["SI"].val),
                        dfdx_MG(dfdx_data["MG"].val)
    {}
    void calculate_G() override {
        phase_free_energy = R*T*(x_CU*log(x_CU) + x_SI*log(x_SI) + x_MG*log(x_MG)) + 2.5*x_CU*x_MG;
    }
    void calculate_dfdx() override {
        dfdx_CU = R*T*(log(x_CU)+1.0)  + 2.5*x_MG;
        dfdx_SI = R*T*(log(x_SI)+1.0);
        dfdx_MG = R*T*(log(x_MG)+1.0)  + 2.5*x_CU;
    }
}

template <int dim, int degree>
class FCC : public PhaseFieldContainer<dim, degree>{
private:
    const scalarValue& x_CU = x_data["CU"].val;
    const scalarValue& x_SI = x_data["SI"].val;
    const scalarValue& x_MG = x_data["MG"].val;

    scalarValue& dfdx_CU = dfdx_data["CU"].val;
    scalarValue& dfdx_SI = dfdx_data["SI"].val;
    scalarValue& dfdx_MG = dfdx_data["MG"].val;
public:
    phase_name = "FCC";
    FCC() : x_CU(x_data["CU"].val),
            x_SI(x_data["SI"].val),
            x_MG(x_data["MG"].val),
            dfdx_CU(dfdx_data["CU"].val),
            dfdx_SI(dfdx_data["SI"].val),
            dfdx_MG(dfdx_data["MG"].val)
    {}
    void calculate_G() override {
        phase_free_energy = // ... ;
    }
    void calculate_dfdx() override {
        dfdx_CU = // ... ;
        dfdx_SI = // ... ;
        dfdx_MG = // ... ;
    }
}
// ...
// Other Phases 
// ...


// 
template <int dim, int degree>
SystemContainer<dim, degree>::SystemContainer(const IsothermalSystem& isoSys, variableContainer<dim,degree,scalarValue>& variable_list){
    uint var_index = 0;
    // For all phase names
    phase_fields.insert({"ExamplePhase", new ExamplePhase<dim,degree>(isoSys.phases.at("ExamplePhase"), variable_list, var_index)});
    phase_fields.insert({"FCC", new FCC<dim,degree>(isoSys.phases.at("FCC"), variable_list, var_index)});
    phase_fields.insert({"BCC", new BCC<dim,degree>(isoSys.phases.at("BCC"), variable_list, var_index)});
    phase_fields.insert({"HCP", new HCP<dim,degree>(isoSys.phases.at("HCP"), variable_list, var_index)});
    // etc.
}

