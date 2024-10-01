// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.
#include "json.hpp"
#include "SystemContainer.h"
void variableAttributeLoader::loadVariableAttributes(){
    #include "IsothermalSystem.h"

    // Open file ------------------------------------------------------------------------
    // Different scope from customPDE, must initialize identical system
    std::string sysFile = "system.json";
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
    // Close the file -------------------------------------------------------------------
    inputFile.close();

    // Make System object
    IsothermalSystem Sys(TCSystem);

    // Load Variable Attributes
    Sys.load_variables(this);
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS([[maybe_unused]] variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 [[maybe_unused]] dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    // --- Get the values and derivatives of the model variables ---
    SystemContainer<dim, degree> sysFields(Sys, variable_list);
    // Solve
    sysFields.initialize_fields();
    sysFields.solve();
    sysFields.submit_fields();
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS([[maybe_unused]] variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 [[maybe_unused]] dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS([[maybe_unused]] variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		[[maybe_unused]] dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}