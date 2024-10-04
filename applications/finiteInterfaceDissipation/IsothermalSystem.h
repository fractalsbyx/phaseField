#include <map>
#include <string>
#include "json.hpp"
#include "Phase.h"
#include "../../include/variableAttributeLoader.h"

struct CompInfo{
    double P; // permeability
};

class IsothermalSystem {
public:
    std::map<std::string, Phase> phases;
    std::map<std::string, CompInfo> comp_info;
    uint N;
    double Vm;
    double eta;
    double Temperature;
    IsothermalSystem(){}
    IsothermalSystem(const nlohmann::json& TCSystem){
        from_json(TCSystem);
    }

    void from_json(const nlohmann::json& j) {
        // Parse Vm and eta
        Vm = j.at("Vm").get<double>();
        eta = j.at("eta").get<double>();

        // Parse components
        for (const auto& [comp_name, comp_data] : j.at("components").items()) {
            CompInfo compInfo;
            compInfo.P = comp_data.at("permeability").get<double>();
            comp_info[comp_name] = compInfo;
        }

        // Parse phases
        for (const auto& [phase_name, phase_data] : j.at("phases").items()) {
            Phase phase;
            phase.name = phase_name;
            phase.sigma = phase_data.at("sigma").get<double>();

            for (const auto& [comp_name, comp_data] : phase_data.items()) {
                if (comp_name != "sigma") {
                    PhaseCompInfo phaseCompInfo;
                    phaseCompInfo.M = comp_data.at("mobility").get<double>();
                    phase.comps[comp_name] = phaseCompInfo;
                }
            }
            phases[phase_name] = phase;
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
            for (const auto& [comp, info] : comp_info){
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
            for (const auto& [comp, info] : comp_info){
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