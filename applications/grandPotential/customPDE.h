#include "../../include/matrixFreePDE.h"
#include <random>

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
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

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
        const unsigned int num_phases = userInputs.get_model_constant_int("num_phases"); // 4
        const unsigned int num_comps  = userInputs.get_model_constant_int("num_comps"); // 3
        const unsigned int num_compVars = num_comps-1;
        
        const std::vector<double> fWell = userInputs.get_model_constant_double_array("fWell");
        const std::vector<std::vector<double>> kWell = reshapeVector(userInputs.get_model_constant_double_array("kWell"),
                                                                num_phases, num_comps);
        const std::vector<std::vector<double>> cmin  = reshapeVector(userInputs.get_model_constant_double_array("cmin"),
                                                                num_phases, num_comps);
        const double L      = userInputs.get_model_constant_double("L");
        const double m0     = userInputs.get_model_constant_double("m0");
        const double kappa  = userInputs.get_model_constant_double("kappa");
        const double gamma  = userInputs.get_model_constant_double("gamma");
        const double M      = userInputs.get_model_constant_double("M");
        const double Va     = userInputs.get_model_constant_double("Va");

        /*
        const std::vector<std::vector<double>> kWell
            {{0.8,0.8,0.8},{1.0,2.0,1.0},{2.0,2.0,2.0},{50.0,50.0,50.0}};
        const std::vector<std::vector<double>> cmin
            {{0.75,0.1,0.15},{0.94,0.03,0.03},{0.67,0.33,0.0},{0.0,0.0,1.0}};
        const std::vector<double> fWell{0.0,0.0,0.0,0.0};
        const double L{1.0}, m0{1.0}, kappa{0.125}, Va{1.0}, gamma{1.5};
        */

        //Declaring random number generator (Type std::mt19937_64)
         engine rng;
         //Declaring distribution (Type std::uniform_real_distribution<double>)
         distribution dist;
	// ================================================================

};
