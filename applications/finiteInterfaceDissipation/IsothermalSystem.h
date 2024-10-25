#ifndef ISOTHERMALSYSTEM_H
#define ISOTHERMALSYSTEM_H

#include "../../include/variableAttributeLoader.h"

#include "Phase.h"
#include "json.hpp"

#include <map> //
#include <string> //
#include <iostream>

struct CompInfo{
    double _P; // permeability
    double P; // permeability
};

class IsothermalSystem {
public:
    std::map<std::string, Phase> phases;
    std::map<std::string, CompInfo> comp_info;
    uint N;
    double length_scale, time_scale, energy_scale;
    double _Vm, Vm;
    double _eta, eta;
    IsothermalSystem(){}
    IsothermalSystem(const nlohmann::json& TCSystem){
        from_json(TCSystem);
    }

    void from_json(const nlohmann::json& j) {
        // Parse Dimensions
        length_scale = j.at("dimensions").at("length_scale").get<double>();
        time_scale = j.at("dimensions").at("time_scale").get<double>();
        energy_scale = j.at("dimensions").at("energy_density_scale").get<double>();

        // Parse Vm and eta
        _Vm = j.at("Vm").get<double>();
        _eta = j.at("eta").get<double>();

        // Parse components
        for (const auto& [comp_name, comp_data] : j.at("components").items()) {
            CompInfo compInfo;
            compInfo._P = comp_data.at("permeability").get<double>();
            comp_info[comp_name] = compInfo;
        }

        // Parse phases
        for (const auto& [phase_name, phase_data] : j.at("phases").items()) {
            Phase phase;
            phase.name = phase_name;
            phase._sigma = phase_data.at("sigma").get<double>();
            phase._mu = phase_data.at("mu").get<double>();

            for (const auto& [comp_name, comp_data] : phase_data.items()) {
                if (comp_name != "sigma" && comp_name != "mu") {
                    PhaseCompInfo phaseCompInfo;
                    phaseCompInfo._M = comp_data.at("mobility").get<double>();
                    phase.comps[comp_name] = phaseCompInfo;
                }
            }
            phases[phase_name] = phase;
        }
        N = phases.size();
        nondimensionalize();
    }

    void nondimensionalize(){
        Vm = _Vm;
        eta = _eta;
        for (auto& [phase_name, phase] : phases){
            phase.mu = phase._mu;
            phase.sigma = phase._sigma;
            for (auto& [comp_name, comp] : phase.comps){
                comp.M = comp._M;
            }
        }
        for (auto& [comp_name, comp] : comp_info){
            comp.P = comp._P;
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
            loader->set_dependencies_gradient_term_RHS(var_index, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
            var_index++;
            for (const auto& [comp, info] : comp_info){
                std::string var_name = phase_name + '_' + comp;
                loader->set_variable_name				(var_index, var_name);
                loader->set_variable_type				(var_index, SCALAR);
                loader->set_variable_equation_type		(var_index, EXPLICIT_TIME_DEPENDENT);
                loader->set_need_value_nucleation		(var_index, true);
                loader->set_dependencies_value_term_RHS (var_index, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
                loader->set_dependencies_gradient_term_RHS(var_index, phase_names+','+grad_phase_names+','+comp_names+','+grad_comp_names);
                var_index++;
            }
        }
        std::cout << "Phase names: " << phase_names << "\n"
                  << "Comp names: " << comp_names << "\n";
    }
};

#endif