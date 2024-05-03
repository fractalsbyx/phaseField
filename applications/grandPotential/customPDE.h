#include "../../include/matrixFreePDE.h"
#include <random>
#include <unordered_map>

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
        //Constructed variables & non-dimensionalization
        divideByX(fWell, Va*m0_dim);
        for(uint i=0;i<kWell.size(); ++i){
            divideByX(kWell[i], Va*m0_dim);
        }
        D = D/(l0*l0);// *m0*Va*Va/(l0*l0);//*m0*l0*l0*l0*l0;
        //Defining seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //Initializing distribution
        dist = distribution(0.0,1.0);
        //Initializing random variable
        rng = engine(seed);
    };

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

    typedef std::mt19937_64 engine;
    typedef std::uniform_real_distribution<double> distribution;


private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif
        #ifdef NUCLEATION_FILE_EXISTS
        double getNucleationProbability(variableValueContainer variable_value, double dV) const;
        #endif
	// Function to set the nucleation probability (in nucleation.h)
  // double getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const

	// ================================================================
	// Methods specific to this subclass
	// ================================================================
    // Function to reshape a 1D vector to a 2D vector with n rows and m columns
    std::vector<std::vector<double>> reshapeVector(const std::vector<double>& inputVector, unsigned int n, unsigned int m) {
        int totalElements = n * m;

        // Check if the total number of elements matches the size of the input vector
        if (inputVector.size() != totalElements) {
            std::cerr << "[customPDE.h] Error: The total number of elements ["
                        << (n*m) << "] does not match the size of the input vector ["
                        << inputVector.size() << "]." << std::endl;
            std::exit(1);
            // Return an empty 2D vector to indicate an error
            return {};
        }

        // Reshape the vector
        std::vector<std::vector<double>> reshapedVector(n, std::vector<double>(m));

        for (int i = 0; i < totalElements; ++i) {
            reshapedVector[i / m][i % m] = inputVector[i];
        }

        return reshapedVector;
    }

    // key is phase ID, value is list of op indeces
    std::unordered_map<unsigned int, std::vector<unsigned int> > make_map(std::vector<int> phase_id){
        std::unordered_map<unsigned int, std::vector<unsigned int>> out;
        for(unsigned int var_index = 0; var_index < phase_id.size(); ++var_index){
            out[phase_id[var_index]].push_back(var_index);
        }
        return out;
    }
    std::unordered_map<unsigned int, unsigned int> make_map_from_name(std::vector<std::string> phase_id, std::unordered_map<std::string, unsigned int> phase_index){
        std::unordered_map<unsigned int, unsigned int> out;
        for(unsigned int var_index = 0; var_index < phase_id.size(); ++var_index){
            out[phase_index[phase_id[var_index]]] = var_index;
        }
        return out;
    }
    std::vector<int> get_phase_index(unsigned int _num_phases, unsigned int _num_ops){
        std::vector<int> out(_num_ops);
        if(boost::iequals(userInputs.get_model_constant_string("set_ids_by"), "INDEX")){
            out = userInputs.get_model_constant_int_array("phase_id");
        }
        else if(boost::iequals(userInputs.get_model_constant_string("set_ids_by"), "NAME")){
            std::unordered_map<std::string, unsigned int> phase_id_by_name;
            std::vector<std::string> phase_names = userInputs.get_model_constant_string_array("phase_names");
            for(unsigned int i = 0; i<_num_phases; ++i){
                phase_id_by_name[phase_names[i]] = i;
            }
            std::vector<std::string> phase_id_names = userInputs.get_model_constant_string_array("phase_id");
            for(unsigned int op = 0; op<_num_ops; ++op){
                out[op] = phase_id_by_name[phase_id_names[op]];
            }
        }
        return out;
    }

    void divideByX(std::vector<double> &vec, double X){
        for(uint i=0; i<vec.size(); ++i){
            vec[i] /= X;
        }
    }

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
        unsigned int num_phases = userInputs.get_model_constant_int("num_phases");
        unsigned int num_comps  = userInputs.get_model_constant_int("num_comps");
        unsigned int num_muFields = num_comps-1;
        unsigned int num_ops = userInputs.get_model_constant_int("num_ops");

        std::vector<double> fWell = userInputs.get_model_constant_double_array("fWell");
        // These variable are read in linearly, but neet to be 2D vectors of shape [num_phases][num_comps-1] 
        std::vector<std::vector<double>> kWell = reshapeVector(userInputs.get_model_constant_double_array("kWell"),
                                                                  num_phases, num_muFields);
        std::vector<std::vector<double>> cmin  = reshapeVector(userInputs.get_model_constant_double_array("cmin"),
                                                                  num_phases, num_muFields);
        double Va     = userInputs.get_model_constant_double("Va");
        double L      = userInputs.get_model_constant_double("L");
        double l0     = userInputs.get_model_constant_double("l0");
        double l_gb   = userInputs.get_model_constant_double("l_gb");
        double sigma  = userInputs.get_model_constant_double("sigma");
        double m0_dim = 6.0*sigma/l_gb;
        //double kappa_dim  = sigma*l_gb*3.0/4.0;
        double m0     = 1.0;
        double kappa  = (1.0/8.0)*(l_gb*l_gb)/(l0*l0); //kappa_dim/(m0_dim*l0*l0);
        double gamma  = userInputs.get_model_constant_double("gamma");
        double D      = userInputs.get_model_constant_double("D");
        
        std::vector<int> phase_index = get_phase_index(num_phases, num_ops);

        // std::unordered_map<unsigned int, std::vector<unsigned int>> ops_of_phase;
        // std::vector<std::string> phase_names = userInputs.get_model_constant_string_array("phase_names");

        //Declaring random number generator (Type std::mt19937_64)
         engine rng;
         //Declaring distribution (Type std::uniform_real_distribution<double>)
         distribution dist;
	// ================================================================

};
