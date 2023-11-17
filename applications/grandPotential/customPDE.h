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


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
        // // Number of componenets 
        // const unsigned int num_comps = 2; // aka K
        // // Number of phases
        // const unsigned int num_phases = 4; // aka N
        // // Constitutional free energy well curvature [phase] [component]
        // const std::vector<std::vector<double>> kWell
        //     {{0.8,0.8,0.8},{1.0,2.0,1.0},{2.0,2.0,2.0},{50.0,50.0,50.0}};
        // const std::vector<double> fWell{0.01,0.0,0.0,0.0};
        // // Composition at constitutional free energy well "minimum"
        // const std::vector<std::vector<double>> cmin
        //     {{0.75,0.1,0.15},{0.94,0.03,0.03},{0.67,0.33,0.0},{0.0,0.0,1.0}};
        // // User input initial (liquid) composition
        // const std::vector<double> c0{0.1, 0.2, 0.7};

        // Number of componenets 
        const unsigned int num_comps = 2; // aka K
        // Number of phases
        const unsigned int num_phases = 4; // aka N
        // Constitutional free energy well curvature [phase] [component]
        const std::vector<std::vector<double>> kWell
            {{0.8,0.8,0.8},{1.0,2.0,1.0},{2.0,2.0,2.0},{50.0,50.0,50.0}};
        const std::vector<double> fWell{0.0,0.0,0.0,0.0};
        // Composition at constitutional free energy well "minimum"
        const std::vector<std::vector<double>> cmin
            {{0.75,0.1,0.15},{0.94,0.03,0.03},{0.67,0.33,0.0},{0.0,0.0,1.0}};
        // User input initial (liquid) composition
        const std::vector<double> c0{0.1, 0.9, 0.7};

        // Initial chemical potential derived from free energy well and c0 assuming a parabolic form
        // mu = k[phase][comp]*(c[comp] - cmin[comp])
        std::vector<std::vector<double>> mu_ini; // set in IC

        const double L{1.0};// mWell{0.1}, kappa{0.0125}, Va{1.0}, gamma{1.5};
        const double M{1.0};

        const double m0 = userInputs.get_model_constant_double("m0"); // energy scale
	    const double gamma = userInputs.get_model_constant_double("gamma"); // interfacial parameter aka alpha
	    const double kappa = userInputs.get_model_constant_double("kappa"); // gradient parameter
        const double Va = userInputs.get_model_constant_double("Va"); // volume scale
        
        //Declaring random number generator (Type std::mt19937_64)
        engine rng;
        //Declaring distribution (Type std::uniform_real_distribution<double>)
        distribution dist;

        
	// ================================================================

};
