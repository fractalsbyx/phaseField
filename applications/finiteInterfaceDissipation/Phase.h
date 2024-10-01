#include <set>
#include <vector>
#include <string>
#include "json.hpp"

struct PhaseCompInfo{
    double M;
};

class Phase {
    public:
    Phase(const nlohmann::json_abi_v3_11_2::detail::iteration_proxy_value<nlohmann::json_abi_v3_11_2::detail::iter_impl<nlohmann::json_abi_v3_11_2::json>> &phase){
        name = phase.key();
        const auto& phase_data = phase.value();

        // Components
        auto comp_info = phase_data["components"];
        std::set<std::string> comps;
        for (auto [name, val] : comp_info){
            double mobility = val.as_object().at("mobility");
            comps.insert({name, {mobility}});
        }
        comps.insert(comp_info.begin(), comp_info.end());
    }
    public:
    std::string name;
    std::map<std::string, PhaseCompInfo>& comps;
};