#include <set>
#include <vector>
#include <string>
#include "json.hpp"

class Phase {
    public:
    Phase(){}
    Phase(const nlohmann::json_abi_v3_11_2::detail::iteration_proxy_value<nlohmann::json_abi_v3_11_2::detail::iter_impl<nlohmann::json_abi_v3_11_2::json>> &phase){
        name = phase.key();
        const auto& phase_data = phase.value();

        // Components
        std::vector<std::string> comp_vec = phase_data["components"];
        std::set<std::string> comps;
        comps.insert(comp_vec.begin(), comp_vec.end());
    }
    public:
    std::string name;
    std::unordered_set<std::string> comps;
};