#include <unordered_map>
#include <string>
#include "json.hpp"
#include "Phase.h"

class IsothermalSystem {
public:
    IsothermalSystem(){}
    IsothermalSystem(nlohmann::json& TCSystem){
        for (const auto& phase : TCSystem["phases"].items()) {
            std::string phase_name = phase.key();
            phases[phase_name] = Phase(phase);
        }
    }
    std::unordered_map<std::string, Phase> phases;
    double Temperature;
};