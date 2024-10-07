// This file will be created or modified by AMMBER
#include "SystemContainer.h"

// Individual phases are derived classes of PhaseFieldContainer
template <int dim, int degree>
class ExamplePhase : public PhaseFieldContainer<dim, degree> {
public:
    const typename PhaseFieldContainer<dim, degree>::scalarValue& x_CU;
    const typename PhaseFieldContainer<dim, degree>::scalarValue& x_SI;
    const typename PhaseFieldContainer<dim, degree>::scalarValue& x_MG;

    typename PhaseFieldContainer<dim, degree>::scalarValue& dfdx_CU;
    typename PhaseFieldContainer<dim, degree>::scalarValue& dfdx_SI;
    typename PhaseFieldContainer<dim, degree>::scalarValue& dfdx_MG;

public:
    ExamplePhase(const IsothermalSystem& isoSys, const std::string& phase_name,
                 const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
                 variableContainer<dim, degree, typename PhaseFieldContainer<dim, degree>::scalarValue>& variable_list) 
        : PhaseFieldContainer<dim, degree>(isoSys, phase_name, phase_fields, variable_list),
          x_CU(this->comp_data.at("CU").x_data.val),
          x_SI(this->comp_data.at("SI").x_data.val),
          x_MG(this->comp_data.at("MG").x_data.val),
          dfdx_CU(this->comp_data.at("CU").dfdx.val),
          dfdx_SI(this->comp_data.at("SI").dfdx.val),
          dfdx_MG(this->comp_data.at("MG").dfdx.val)
    {}

    void calculate_G() override {
        this->phase_free_energy = (x_CU * log(x_CU) + x_SI * log(x_SI) + x_MG * log(x_MG)) + 2.5 * x_CU * x_MG;
    }

    void calculate_dfdx() override {
        dfdx_CU = (log(x_CU) + 1.0) + 2.5 * x_MG;
        dfdx_SI = (log(x_SI) + 1.0);
        dfdx_MG = (log(x_MG) + 1.0) + 2.5 * x_CU;
    }
};



//template <int dim, int degree>
//class FCC : public PhaseFieldContainer<dim, degree>{
//private:
//    const scalarValue& x_CU;
//    const scalarValue& x_SI;
//    const scalarValue& x_MG;
//
//    scalarValue& dfdx_CU;
//    scalarValue& dfdx_SI;
//    scalarValue& dfdx_MG;
//    scalarValue& d2fdx2_CU;
//    scalarValue& d2fdx2_SI;
//    scalarValue& d2fdx2_MG;
//public:
//    phase_name = "FCC";
//    FCC(const IsothermalSystem& isoSys, const std::string& phase_name,
//        const std::map<std::string, PhaseFieldContainer<dim, degree>*>& phase_fields,
//        variableContainer<dim,degree,scalarValue>& variable_list) : 
//        x_CU(x_data["CU"].val),
//        x_SI(x_data["SI"].val),
//        x_MG(x_data["MG"].val),
//        dfdx_CU(dfdx_data["CU"].val),
//        dfdx_SI(dfdx_data["SI"].val),
//        dfdx_MG(dfdx_data["MG"].val)
//    {}
//    void calculate_G() override {
//        phase_free_energy = // ... ;
//    }
//    void calculate_dfdx() override {
//        dfdx_CU = // ... ;
//        dfdx_SI = // ... ;
//        dfdx_MG = // ... ;
//    }
//}
// ...
// Other Phases 
// ...

// 
template <int dim, int degree>
SystemContainer<dim, degree>::SystemContainer(const IsothermalSystem& _isoSys, variableContainer<dim,degree,scalarValue>& _variable_list) :
        isoSys(_isoSys), variable_list(_variable_list){
    // For all phase names
    phase_fields.insert({"ExamplePhase", new ExamplePhase<dim,degree>(isoSys, "ExamplePhase", phase_fields, variable_list)});
    //phase_fields.insert({"FCC", new FCC<dim,degree>(isoSys, "FCC", phase_fields, variable_list)});
    //phase_fields.insert({"BCC", new BCC<dim,degree>(isoSys, "BCC", phase_fields, variable_list)});
    //phase_fields.insert({"HCP", new HCP<dim,degree>(isoSys, "HCP", phase_fields, variable_list)});
    // etc.
}
