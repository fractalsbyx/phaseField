#ifndef PHASE_H
#define PHASE_H


#include <set>
#include <vector>
#include <string>

struct PhaseCompInfo{
    double M;
};

struct Phase {
public:
    Phase(){}
    /*
    Phase(const nlohmann::json &phases, const std::string& phase_name){
        name = phase_name;
        const auto& phase_data = phases[phase_name];

        // Constants
        sigma = phase_data["sigma"];

        // Components
        nlohmann::json comp_info = phase_data["components"];
        std::set<std::string> comps;
        for (auto [comp_name, val] : comp_info){
            double mobility = comp_info[comp_name]["mobility"];
            comps.insert({comp_name, {mobility}});
        }
        comps.insert(comp_info.begin(), comp_info.end());
    }
    */
    std::string name;
    std::map<std::string, PhaseCompInfo> comps;
    double sigma;
    double mu;
};

#endif