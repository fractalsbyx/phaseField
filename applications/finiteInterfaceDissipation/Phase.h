#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <list>
#include <string>
#include <boost/algorithm/string.hpp>
#include "json.hpp"

struct InteractionParameter {
    std::vector<std::vector<std::string>> occupation;
    uint degree;
    double value;
};

class Phase {
    public:
    Phase(){}
    Phase(const nlohmann::json_abi_v3_11_2::detail::iteration_proxy_value<nlohmann::json_abi_v3_11_2::detail::iter_impl<nlohmann::json_abi_v3_11_2::json>> &phase){
        name = phase.key();
        const auto& phase_data = phase.value();

        // Number of sublattices
        num_sublattices = phase_data["num_sublattices"];

        // Number of Sites
        num_sites = phase_data["num_sites"].get<std::vector<double>>();

        // Is interstitial
        is_interstitial = phase_data["num_sites"].get<std::vector<bool>>();

        // Total number of sites
        total_num_sites = 0;
        total_num_subst_sites = 0;
        for (uint s = 0; s < num_sublattices; ++s){
            total_num_sites += num_sites[s];
            total_num_subst_sites += num_sites[s]*(!is_interstitial[s]);
        }

        // Occupation (sublattice_comps)
        std::vector<std::unordered_set<std::string>> sublattice_vec;
        for (const auto& sublattice : phase_data["occupation"]) {
            std::unordered_set<std::string> sublattice_set(sublattice.begin(), sublattice.end());
            sublattice_vec.push_back(sublattice_set);
            comps.insert(sublattice_set.begin(), sublattice_set.end());
        }
        sublattice_comps = sublattice_vec;

        // Interaction parameters
        for (const auto& param : phase_data["interaction parameters"]) { //.items() for key-val pair
            InteractionParameter L;
            for (const auto& sublattice : param["occupation"]) {
                std::vector<std::string> constituents(sublattice.begin(), sublattice.end());
                L.occupation.push_back(constituents);
            }
            L.degree = param["degree"];
            L.value = param["value"];
        }
    }
    public:
    std::string name;
    uint num_sublattices;
    std::vector<double> num_sites;
    std::vector<double> a;
    std::vector<bool> is_interstitial;
    std::list<uint> mixed_sublattices;
    double total_num_sites;
    double total_num_subst_sites;
    std::vector<std::set<std::string>> sublattice_comps;
    std::unordered_set<std::string> comps;
    std::unordered_map<std::string, InteractionParameter> L_parameter;
};