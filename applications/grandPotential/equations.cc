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
//#include "customPDE.h"
void variableAttributeLoader::loadVariableAttributes(){
    const unsigned int num_phases = 3;
    const unsigned int num_comps = 2;
    std::string list_valn = "";
    std::string list_gradn = "";
    std::string list_valmu = "";
    std::string list_gradmu = "";
    // std::string list_valdndt = "";
    for (unsigned int var_index=0; var_index<num_phases; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));
        list_valn.append(var_name+", ");
        list_gradn.append("grad("+var_name+"), ");
        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,EXPLICIT_TIME_DEPENDENT);
    }
    for (unsigned int comp_index=0; comp_index<num_comps; comp_index++){
        unsigned int var_index = num_phases+comp_index;
        std::string var_name = "mu";
        var_name.append(std::to_string(comp_index));
        list_valmu.append(var_name+", ");
        list_gradmu.append("grad("+var_name+"), ");
        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,EXPLICIT_TIME_DEPENDENT);
    }

    std::cout << list_valn << "| " << list_valmu << "| "
        << list_gradn << "| " << list_gradmu << "\n";
    std::string dep_valn = list_valn+list_gradn+list_valmu+list_gradmu;
    std::string dep_valmudndt = list_valn+list_gradn+list_valmu+list_gradmu;
    dep_valn.pop_back();        // remove last comma and space
    dep_valmudndt.pop_back();
    dep_valn.pop_back();
    dep_valmudndt.pop_back();
    for (unsigned int var_index=0; var_index<num_phases; var_index++){ // phase order params
        set_dependencies_value_term_RHS(var_index, dep_valn); // order params, saved dndt vals, mu vals
        set_dependencies_gradient_term_RHS(var_index, "");
    }
    for (unsigned int var_index=num_phases; var_index<num_phases+num_comps; var_index++){ // mu vals
        set_dependencies_value_term_RHS(var_index, dep_valmudndt);
        set_dependencies_gradient_term_RHS(var_index, list_gradmu);
    }
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
//std::cout << "Explicit";
std::vector<scalarvalueType> eta_values(num_phases);
std::vector<scalargradType> eta_gradients(num_phases);
std::vector<scalarvalueType> dndt_values(num_phases);
std::vector<scalargradType> dndt_grad(num_phases);
std::vector<scalarvalueType> mu_values(num_comps);
std::vector<scalargradType> mu_gradients(num_comps);
scalarvalueType sum_nsq = constV(0.0);
std::vector<scalarvalueType> omegaC(num_phases);
// get fields
for (unsigned int i=0; i<num_phases; ++i){
	eta_values[i] = variable_list.get_scalar_value(i);
    eta_gradients[i] = variable_list.get_scalar_gradient(i);
}
for (unsigned int i=0; i<num_comps; ++i){
	mu_values[i] = variable_list.get_scalar_value(i+num_phases);
	mu_gradients[i] = variable_list.get_scalar_gradient(i+num_phases);
}
// calc sum_nsq
for (unsigned int i=0; i<num_phases; ++i){
	sum_nsq += eta_values[i]*eta_values[i];
}
// START COPIED PORTION
for (unsigned int i=0; i<num_phases; ++i){
    omegaC[i] = constV(fWell[i]);
    for (unsigned int j=0; j<num_comps; ++j){
        omegaC[i] -= constV(0.5)*mu_values[j]*mu_values[j]/constV(Va*Va*kWell[i][j])
            + mu_values[j]*constV(cmin[i][j]/Va);
    }
}
// calc weak form for eta
for (unsigned int i=0; i < num_phases; ++i){
    dndt_values[i] += m0*(eta_values[i]*eta_values[i]*eta_values[i] - eta_values[i]);
    dndt_grad[i] += kappa*eta_gradients[i];
    for (unsigned int j=0; j<num_phases; ++j){
        if(i != j){ // interface term
            dndt_values[i] += constV(m0*2.0)*eta_values[i]*gamma*eta_values[j]*eta_values[j];
        } else { // dhdn stuff
            dndt_values[i] += constV(2.0)*omegaC[i]*eta_values[i]/sum_nsq;
        }
        dndt_values[i] -= constV(2.0)*omegaC[j]*eta_values[j]*(eta_values[j]*eta_values[i])/(sum_nsq*sum_nsq); // dhdn stuff
    }
}
// END COPIED PORTION

// calc h and dhdn
std::vector<scalarvalueType> h(num_phases);
std::vector<std::vector<scalarvalueType>> dhdn(num_phases, std::vector<scalarvalueType> (num_phases, constV(0.0)));
//std::vector<scalarvalueType> sum_dhdn_omega(num_phases);
for (unsigned int i=0; i<num_phases; ++i){
    scalarvalueType eta_i_sq = eta_values[i]*eta_values[i];
    h[i] = eta_i_sq/sum_nsq;
    for (unsigned int j=0; j<num_phases; ++j){
        if(i == j){
            dhdn[i][j] += 2.0*eta_values[i]/(sum_nsq);
            //sum_dhdn_omega[i] += omegaC[i]*2*eta_values[i]*(sum_nsq-eta_i_sq);
        }
        else{
            //sum_dhdn_omega[i] -= omegaC[j]*2*eta_values[j]*(eta_i_sq);
        }
        dhdn[i][j] -= 2.0*eta_values[i]*eta_values[i]*eta_values[j]/(sum_nsq*sum_nsq);
    }
    //sum_dhdn_omega[i] /= sum_nsq*sum_nsq;
}

// calc weak form for mu
std::vector<scalarvalueType> dmudtValue(num_comps);
std::vector<scalargradType> dmudtGrad(num_comps);
//std::vector<scalarvalueType> chiAA(num_comps, constV(0.0));
for (unsigned int i=0; i<num_comps; ++i){
    scalarvalueType susceptibility = constV(0.0);
    for (unsigned int j=0; j<num_phases; ++j){
        susceptibility += h[j]/(Va*Va*kWell[j][i]);
        //chiAA[i] += (eta_values[j]*eta_values[j]/sum_nsq)/(Va*Va*kWell[j][i]);
    }
    dmudtGrad[i] -= M*mu_gradients[i]/susceptibility;
    // xan // dmudtGrad[i] = M*mu_gradients[i]/chiAA[i];
    for (unsigned int j=0; j<num_phases; ++j){
        for (unsigned int k=0; k<num_phases; ++k){
            scalarvalueType drhodn = dhdn[k][j]*(mu_values[i]/(Va*kWell[k][i]) + constV(cmin[k][i]));
            dmudtValue[i] -= dndt_values[j]*drhodn;
        }
    }
}

// submit terms
for (unsigned int i=0; i<num_phases; ++i){
    variable_list.set_scalar_value_term_RHS(i,eta_values[i] +
        dndt_values[i]*constV(userInputs.dtValue));
    variable_list.set_scalar_gradient_term_RHS(i,dndt_grad[i]*constV(userInputs.dtValue));
}
for (unsigned int i=0; i<num_comps; ++i){
    dmudtValue[i] = mu_values[i] + dmudtValue[i]*constV(userInputs.dtValue);
    dmudtGrad[i] = dmudtGrad[i]*constV(userInputs.dtValue);
    variable_list.set_scalar_value_term_RHS(i+num_phases,dmudtValue[i]);
    variable_list.set_scalar_gradient_term_RHS(i+num_phases,dmudtGrad[i]);
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
