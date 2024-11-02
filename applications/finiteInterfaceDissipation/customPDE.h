#include "../../include/matrixFreePDE.h"
#include <random>
//#include "customPhases.cc"
#include "SystemContainer.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {
        std::cout << "Starting customPDE initializer...\n";
        //Read system json file to TCSystem
        parseSystem();
        // Zero grad constant
        ZERO = 0.0*ZERO;
        //Defining seed
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //Initializing distribution
        dist = distribution(0.0,1.0);
        //Initializing random variable
        rng = engine(seed);
        std::cout << "customPDE initialized\n";
    };

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC) override;

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC) override;

    typedef std::mt19937_64 engine;
    typedef std::uniform_real_distribution<double> distribution;


private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;
	#endif
    // Function to set the nucleation probability (in nucleation.h)
    #ifdef NUCLEATION_FILE_EXISTS
    double getNucleationProbability(variableValueContainer variable_value, double dV) const;
    double getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const;
    #endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================

    void parseSystem(){
        // Create an input file stream
        std::ifstream inputFile(sysFile);
        // Check if the file was successfully opened
        if (!inputFile.is_open()) {
            std::cerr << "Could not open the file: " << sysFile << std::endl;
            std::exit(1);
        }
        // Parse the JSON file into a JSON object
        nlohmann::json TCSystem;
        try {
            inputFile >> TCSystem;
        } catch (nlohmann::json::parse_error& e) {
            std::cerr << "JSON parse error: " << e.what() << std::endl;
            std::exit(1);
        }
        // Close the file
        inputFile.close();

        Sys = IsothermalSystem(TCSystem);
    }

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
        IsothermalSystem Sys;
        // File for system name
        std::string sysFile = "system.json";// userInputs.get_model_constant_string("sysFile");

        double T;
        double sigma;
        double l_gb;

        // prm constants
        double r0 = userInputs.get_model_constant_double("r0");

        // Zero vector
        scalargradType ZERO;
        // pi
        const double pi = 3.141592653589793238;
        // Ideal gas
        const double R = 8.314;
        //Declaring random number generator (Type std::mt19937_64)
        engine rng;
        //Declaring distribution (Type std::uniform_real_distribution<double>)
        distribution dist;
	// ================================================================

};