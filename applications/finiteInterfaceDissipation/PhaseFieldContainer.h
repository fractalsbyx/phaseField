#include "Phase.h"
#include <deal.II/base/vectorization.h>
//#include "../../include/variableContainer.h"


template <int dim>
struct FieldContainer{
    dealii::VectorizedArray<double> val;
    dealii::Tensor<1, dim, dealii::VectorizedArray<double>> grad;
};
struct InteractionData{
        dealii::VectorizedArray<double> full_product;
        std::vector<dealii::VectorizedArray<double>> excluded_product;
};
template <int dim, int degree>
class PhaseFieldContainer{
    typedef dealii::VectorizedArray<double> scalarValue;
    typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
    #define constV(a) dealii::make_vectorized_array(a)

    PhaseFieldContainer(const Phase& phase, variableContainer<dim,degree,scalarValue>& variable_list, uint& var_index) :
                                    name(phase.name),
                                    num_sublattices(phase.num_sublattices),
                                    num_sites(phase.num_sites),
                                    a(phase.a),
                                    is_interstitial(phase.is_interstitial),
                                    mixed_sublattices(phase.mixed_sublattices),
                                    total_num_sites(phase.total_num_sites),
                                    total_num_subst_sites(phase.total_num_subst_sites),
                                    sublattice_comps(phase.sublattice_comps),
                                    comps(phase.comps),
                                    L_parameter(phase.L_parameter),
                                    variable_list(variable_list),
                                    first_var_index(var_index){
        y_data = std::vector<std::unordered_map<std::string, FieldContainer<dim>>>  (num_sublattices);
        dFdy_data = std::vector<std::unordered_map<std::string, FieldContainer<dim>>>  (num_sublattices);
        dydt_data = std::vector<std::unordered_map<std::string, FieldContainer<dim>>>  (num_sublattices);
        for (unsigned int s=0; s < num_sublattices; s++){
            double a_s = num_sites[s] / total_num_sites;
            if (sublattice_comps[s].size() > 1){ // saving memory for fixed sublattices
            for (const std::string& comp : phase.sublattice_comps[s]){
                y_data[s][comp].val = variable_list.get_scalar_value(var_index);
                y_data[s][comp].grad = variable_list.get_scalar_gradient(var_index++);
                c_phase[comp].val += a_s * y_data[s][comp].val;
                c_phase[comp].grad += a_s * y_data[s][comp].grad;
            }
            }
            else{
                auto comp = *(phase.sublattice_comps[s].begin());
                y_data[s][comp].val = constV(1.0);
                c_phase[comp].val += constV(a_s);
            }
        }
        phi_phase.val *= 0.0;
        phi_phase.grad *= 0.0;
    }

    ~PhaseFieldContainer(){}

    public:

    scalarValue dimensionless_entropy() {
        scalarValue S_phase = constV(0.0);
        for (uint s = 0; s < phase.sublattice_comps.size(); ++s){
            double a = phase.num_sites[s]/phase.total_num_sites;
            for (const std::string& constituent : phase.sublattice_comps[s]){
                S_phase += a * y_data[s][constituent] * std::log(y_data[s][constituent]);
            }
        }
        return S_phase;
    }

    void calcEntropicMu() {
        for (uint s = 0; s < sublattice_comps.size(); ++s){
            for (const std::string& constituent : sublattice_comps[s]){
                dFdy_data[s][constituent].val += a[s] * (constV(1.0) + std::log(y_data[s][constituent].val));
                dFdy_data[s][constituent].grad += a[s] * y_data[s][constituent].grad/y_data[s][constituent].val;
            }
        }
    }

    scalarValue SublatticeTerm(std::vector<std::string>& constituents, uint sublattice){
        switch (constituents.size()){
            case 1:
                // Single element
                return (constituents[0] == "*") ? constV(1.0) : y_data[sublattice][constituents[0]].val;
            break;

            case 2:
                // Redlich-Kister
                return y_data[sublattice][constituents[0]].val - y_data[sublattice][constituents[1]].val;
            break;
            
            default:
                // TODO: Muggianu
                return constV(1.0);
        }
    }

    InteractionData Interaction(const InteractionParameter& L){
        const auto& occupation = L.occupation;
        std::vector<scalarValue> terms          (num_sublattices, constV(1.0));
        std::vector<scalarValue> othersProduct  (num_sublattices, constV(1.0));
        scalarValue prod_term = constV(1.0);
        for (uint s = 0; s < num_sublattices; ++s){
            scalarValue sl_term = SublatticeTerm(occupation[s], s);
            prod_term *= sl_term;
            for(uint s1 = 0; s1 < num_sublattices; ++s1){
                othersProduct[s1] *= (s1!=s) ? sl_term : 1.0;
            }
        }
        return {prod_term, othersProduct};
    }

    void updateMu(const InteractionParameter& L, const InteractionData& id) {
        const auto& occupation = L.occupation;
        for (uint s = 0; s < occupation.size(); ++s){
            auto& constituents = occupation[s];
            scalarValue diff;
            switch (constituents.size()){
                case 1:
                    // Single element
                    if (constituents[0] == "*"){
                        dFdy_data[s][constituents[0]].val = id.excluded_product[s];
                    }
                break;

                case 2:
                    // Redlich-Kister
                    diff = (y_data[s][constituents[0]].val - y_data[s][constituents[1]].val);
                    dFdy_data[s][constituents[0]].val += id.excluded_product[s] * (double)L.degree * dealii::Utilities::pow(diff, L.degree-1);
                    dFdy_data[s][constituents[1]].val +=-id.excluded_product[s] * (double)L.degree * dealii::Utilities::pow(diff, L.degree-1);
                    dFdy_data[s][constituents[0]].grad += id.excluded_product[s] * (double)(L.degree * (L.degree-1)) * dealii::Utilities::pow(diff, L.degree-2) * y_data[s][constituents[0].grad];
                    dFdy_data[s][constituents[1]].grad += id.excluded_product[s] * (double)(L.degree * (L.degree-1)) * dealii::Utilities::pow(diff, L.degree-2) * y_data[s][constituents[1].grad];
                break;

                default:
                    // TODO: Muggianu
                    dFdy_data[s][constituents[0]].val *= 1.0;
            }
        }
    }

    void calculateMuBar(){
        for(const auto& comp : comps){
            mu_bar[comp].val = constV(0.0);
            mu_bar[comp].grad *= constV(0.0);
        }
        for(uint s = 0; s < num_sublattices; ++s){
            for(const std::string& constituent : sublattice_comps[s]){
                mu_bar[constituent].val +=  a[s] * dFdy_data[s][constituent].val;
                mu_bar[constituent].grad += a[s] * dFdy_data[s][constituent].grad;
            }
        }
    }

    void calculateFlux(){
        scalarGrad temp;
        temp *= 0.0;
        const scalarValue VAL_ZERO = constV(0.0);
        const scalarGrad VEC_ZERO = temp;


        scalarValue denom_sum;
        std::unordered_map<std::string, FieldContainer> sum_quotient;
        FieldContainer forward_flux_sum;
        FieldContainer backflux_term;

        for(const auto& comp : comps){
            denom_sum = VAL_ZERO;
            forward_flux_sum.val = VAL_ZERO;
            forward_flux_sum.grad = VEC_ZERO;
            for(const auto& s : sublattices_containing[comp]){
                dydt_data[s][comp].grad     = K_flux[s][comp] * (-dFdy_data[s][comp].grad + mu_bar[comp].grad);
                forward_flux_sum.val  += K_flux[s][comp] * (-dFdy_data[s][comp].val  + mu_bar[comp].val ); //complete after this loop
                forward_flux_sum.grad += K_flux[s][comp] * (-dFdy_data[s][comp].grad + mu_bar[comp].grad); //complete after this loop
                for(const std::string& constituent : sublattice_comps[s]){
                    denom_sum += a[s] * K_flux[s][constituent] * y_data[s][constituent].val; // complete after s loop
                }
            }
            sum_quotient[comp].val = forward_flux_sum.val / denom_sum;
            sum_quotient[comp].grad= forward_flux_sum.grad / denom_sum;
        }
        for(uint s = 0; s < num_sublattices; ++s){
            for(const std::string& constituent : sublattice_comps[s]){
                backflux_term.val = VAL_ZERO;
                backflux_term.grad = VEC_ZERO;
                for(const std::string& comp : sublattice_comps[s]){
                    backflux_term.val += sum_quotient[comp].val;
                    backflux_term.grad+= sum_quotient[comp].grad;
                }
                backflux_term.val   *= a[s] * y_data[s][constituent].val * K_flux[s][constituent];
                backflux_term.grad  *= a[s] * y_data[s][constituent].val * K_flux[s][constituent];
                dydt_data[s][constituent].val  -= backflux_term.val;
                dydt_data[s][constituent].grad -= backflux_term.grad;
            }
        }
    }

    void submit_dydt(double dt){
        uint var_index = first_var_index;
        for (const uint& s : mixed_sublattices){
            for (const std::string& constituent : phase.sublattice_comps[s]){
                variable_list.set_scalar_value_term_RHS(var_index,   y_data[s][constituent].val +    dt * dydt_data[s][constituent].val);
                variable_list.set_scalar_value_term_RHS(var_index++,                                -dt * dydt_data[s][constituent].grad);
            }
        }
    }


    protected:
    //references to phase object
    std::string& name;
    uint& num_sublattices;
    std::vector<double>& num_sites;
    std::vector<double>& a;
    std::vector<bool>& is_interstitial;
    std::list<uint>& mixed_sublattices;
    double& total_num_sites;
    double& total_num_subst_sites;
    std::vector<std::set<std::string>>& sublattice_comps;
    std::unordered_set<std::string>& comps;
    std::unordered_map<std::string, InteractionParameter>& L_parameter;
    variableContainer<dim,degree,scalarValue>& variable_list;
    uint first_var_index;

    std::unordered_map<std::string, uint>& sublattices_containing;

    std::unordered_map<std::string, double>& K_flux;
    std::unordered_map<std::string, double>& K_sub;

    std::vector<std::unordered_map<std::string, FieldContainer<dim>>> y_data;
    std::vector<std::unordered_map<std::string, FieldContainer<dim>>> dFdy_data;
    std::vector<std::unordered_map<std::string, FieldContainer<dim>>> dydt_data;

    std::unordered_map<std::string, FieldContainer<dim>> mu_bar;
    std::unordered_map<std::string, FieldContainer<dim>> c_phase;
    FieldContainer<dim> phi_phase;
    scalarValue G;
};
