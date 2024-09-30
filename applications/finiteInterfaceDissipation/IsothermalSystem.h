#include <map>
#include <string>
#include "json.hpp"
#include "Phase.h"
#include "../../include/variableAttributeLoader.h"

struct CompInfo{
    double P; // permeability
};

class IsothermalSystem {
    std::map<std::string, Phase> phases;
    std::map<std::string, CompInfo> comps;
    double Temperature;
public:
    IsothermalSystem(){}
    IsothermalSystem(nlohmann::json& TCSystem){
        for (const auto& phase : TCSystem["phases"].items()) {
            std::string phase_name = phase.key();
            phases[phase_name] = Phase(phase);
        }
        for (const auto& comp : TCSystem["components"].items()) {
            CompInfo comp_info;
            comp_info.P = comp.value()["permeability"];
            comps.insert({comp.key(), comp_info});
        }
    }

    void load_variables(variableAttributeLoader* loader){
        // Get names for composition fields
        std::string phase_names;
        std::string grad_phase_names;
        std::string comp_names;
        std::string grad_comp_names;
        for (const auto& [phase_name, phase] : phases){
            phase_names.append(phase_name + ",");
            grad_phase_names.append("grad(" + phase_name + "),");
            for (const std::string& comp : comps){
                std::string var_name = phase_name + '_' + comp;
                comp_names.append(var_name + ",");
                grad_comp_names.append("grad("+var_name+"),");
            }
        }
        // remove commas
        phase_names.pop_back(); 
        grad_phase_names.pop_back();
        comp_names.pop_back();
        grad_comp_names.pop_back();

        uint var_index = 0;
        // Assign fields
        for (const auto& [phase_name, phase] : phases){
            loader->set_variable_name				(var_index, phase_name);
            loader->set_variable_type				(var_index, SCALAR);
            loader->set_variable_equation_type		(var_index, EXPLICIT_TIME_DEPENDENT);
            // loader->set_allowed_to_nucleate			(var_index, (phase_index>0));
            loader->set_need_value_nucleation		(var_index, true);
            loader->set_dependencies_value_term_RHS (var_index, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
            loader->set_dependencies_gradient_term_RHS(var_index++, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
            for (const std::string& comp : comps){
                std::string var_name = phase_name + '_' + comp;
                loader->set_variable_name				(var_index, var_name);
                loader->set_variable_type				(var_index, SCALAR);
                loader->set_variable_equation_type		(var_index, EXPLICIT_TIME_DEPENDENT);
                loader->set_need_value_nucleation		(var_index, true);
                loader->set_dependencies_value_term_RHS (var_index, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
                loader->set_dependencies_gradient_term_RHS(var_index++, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
            }
        }
        std::cout << "Finished loadVariableAttributes\n";
    }

};