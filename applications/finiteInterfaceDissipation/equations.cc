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
#include "PhaseFieldContainer.h"
void variableAttributeLoader::loadVariableAttributes(){
    #include "IsothermalSystem.h"
    uint num_ops = 4;

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
    // Close the file
    inputFile.close();
    IsothermalSystem Sys(TCSystem);

    // Get names for composition fields
    std::string comp_names;
    std::string grad_comp_names;
    for (const auto& [phase_name, phase] : Sys.phases){
        for (unsigned int s=0; s<phase.num_sublattices; s++){
            if (phase.sublattice_comps[s].size()>1){ // saving memory for fixed sublattices
                for (const std::string& comp : phase.sublattice_comps[s]){
                    std::string var_name = phase_name + '_' + std::to_string(s) + '_' + comp;
                    comp_names.append(var_name);
                    grad_comp_names.append("grad("+var_name+")");
                }
            }
        }
    }
    comp_names.pop_back(); // remove comma
    grad_comp_names.pop_back();

    // Get names for order parameter fields
    std::string op_names;
    std::string grad_op_names;
    for (unsigned int i=0; i<num_ops; ++i){
        std::string var_name = "phi" + std::to_string(i);
        op_names.append(var_name);
        grad_op_names.append("grad("+var_name+")");
    }
    op_names.pop_back(); // remove comma
    grad_op_names.pop_back();

    uint var_index = 0;
    // Assign composition fields
    for (const auto& [phase_name, phase] : Sys.phases){
        for (unsigned int s=0; s<phase.num_sublattices; s++){
            if (phase.sublattice_comps[s].size()>1){ // saving memory for fixed sublattices
                for (const std::string& comp : phase.sublattice_comps[s]){
                    std::string var_name = phase_name + '_' + std::to_string(s) + '_' + comp;
                    set_variable_name				(var_index, var_name);
                    set_variable_type				(var_index, SCALAR);
                    set_variable_equation_type		(var_index, EXPLICIT_TIME_DEPENDENT);
                    set_need_value_nucleation		(var_index, true);
                    set_dependencies_value_term_RHS (var_index, op_names+','+grad_op_names+','+comp_names+','+grad_comp_names);
                    set_dependencies_gradient_term_RHS(var_index++, op_names+','+grad_op_names+','+comp_names+','+grad_comp_names);
                }
            }
        }
    }

    // Assign order parameter fields
    for (unsigned int i=0; i<num_ops; ++i){
        std::string var_name = "phi" + std::to_string(i);
        set_variable_name				(var_index, var_name);
        set_variable_type				(var_index, SCALAR);
        set_variable_equation_type		(var_index, EXPLICIT_TIME_DEPENDENT);
        // set_allowed_to_nucleate			(var_index, (op_index>0));
        set_need_value_nucleation		(var_index, true);
        set_dependencies_value_term_RHS (var_index, op_names+','+grad_op_names+','+comp_names+','+grad_comp_names);
        set_dependencies_gradient_term_RHS(var_index++, op_names+','+grad_op_names+','+comp_names+','+grad_comp_names);
    }
    std::cout << "Finished loadVariableAttributes\n";
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
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---
// Order Parameters
std::vector<scalarvalueType> phi_val(num_ops);
std::vector<scalargradType> phi_grad(num_ops);
// 
std::vector<scalarvalueType> dFdphi_val(num_ops, constV(0.0));
std::vector<scalargradType> dFdphi_grad(num_ops);
// 
std::vector<scalarvalueType> dphidt_val(num_ops, constV(0.0));
std::vector<scalargradType> dphidt_grad(num_ops);
// Local composition [component_name]
std::unordered_map<std::string, scalarvalueType> c_val;
std::unordered_map<std::string, scalargradType> c_grad;

PhaseFieldContainer<dim, degree> phases()
// Retrieve fields
unsigned int var_index = 0;
for (const auto& phase : phases){
    phase.
}
for (unsigned int i=0; i<num_ops; ++i){
	phi_val[i] =  variable_list.get_scalar_value(var_index);
	phi_grad[i] = variable_list.get_scalar_gradient(var_index++);
    std::string phase = op_phase_name[i];
    phi_phase_val[phase] += phi_val[i];   // This and phase comps need order parameters to sum to 1,
    phi_phase_grad[phase] += phi_grad[i]; // otherwise this must be done later with interpolation function.
    for (const std::string& comp : Sys.phases.at(phase).comps){
        c_val[comp] += c_phase_val[phase][comp] * phi_val[i];
        c_grad[comp] +=  c_phase_grad[phase][comp] * phi_val[i] + c_phase_val[phase][comp] * phi_grad[i];
    }
}

// Calculate G for each phase
for (const auto& [key, phase] : Sys.phases){
    scalarvalueType& G_phase = G[phase.name];
    //scalarvalueType 
    G_phase = constV(0.0);
    G_phase += Sys.Temperature*entropy(phase, y_val[phase.name]);
    calcEntropicMu(phase, y_val[phase.name], y_grad[phase.name], mu_val[phase.name], mu_grad[phase.name]);
    for (auto& [key2, parameter] : phase.L_parameter){
        InteractionData L_data = Interaction(parameter, y_val[phase.name]);
        G_phase += parameter.value * L_data.full_product;
        updateMu(parameter, L_data, y_val[phase.name], y_grad[phase.name], mu_val[phase.name], mu_grad[phase.name]);
    }
}

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
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
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
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}