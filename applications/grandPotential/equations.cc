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

void variableAttributeLoader::loadVariableAttributes(){
    const unsigned int num_ops{3};
    const unsigned int num_muFields{2};
    std::string string_valn = "";
    std::string string_valdndt = "";
    std::string string_valmu = "";
    std::string string_gradn = "";
    std::string string_gradmu = "";
for (unsigned int var_index=0; var_index<num_ops; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));
        string_valn.append(var_name+",");
        string_gradn.append("grad("+var_name+"),");
        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,EXPLICIT_TIME_DEPENDENT);
        
    }
    for (unsigned int var_index=0; var_index<num_muFields; var_index++){
        std::string var_name = "mu";
        var_name.append(std::to_string(var_index));
        string_valmu.append(var_name+",");
        string_gradmu.append("grad("+var_name+"),");
        set_variable_name				(num_ops+var_index,var_name);
    	set_variable_type				(num_ops+var_index,SCALAR);
    	set_variable_equation_type		(num_ops+var_index,EXPLICIT_TIME_DEPENDENT);
    }
    for (unsigned int var_index=0; var_index<num_ops; var_index++){
        std::string var_name = "dndt";
        var_name.append(std::to_string(var_index));
        string_valdndt.append(var_name+",");
        set_variable_name				(num_ops+num_muFields+var_index,var_name);
    	set_variable_type				(num_ops+num_muFields+var_index,SCALAR);
    	set_variable_equation_type		(num_ops+num_muFields+var_index,AUXILIARY);
    }
    std::cout << string_valn << " | " << string_valmu << " | " << string_valdndt << " | "
        << string_gradn << " | " << string_gradmu << "\n";
    std::string dep_valn = string_valn+string_valdndt;
    std::string dep_valmudndt = string_valn+string_valmu+string_valdndt;
    dep_valn.pop_back();
    dep_valmudndt.pop_back();
    string_gradmu.pop_back();
    string_gradn.pop_back();
    for (unsigned int var_index=0; var_index<num_ops; var_index++){
        set_dependencies_value_term_RHS(var_index, dep_valn);
        set_dependencies_gradient_term_RHS(var_index, "");
    }
    for (unsigned int var_index=num_ops; var_index<num_ops+num_muFields; var_index++){
        set_dependencies_value_term_RHS(var_index, dep_valmudndt);
        set_dependencies_gradient_term_RHS(var_index, string_gradmu);
    }
    for (unsigned int var_index=num_ops+num_muFields; var_index<2*num_ops+num_muFields; var_index++){
        set_dependencies_value_term_RHS(var_index, dep_valmudndt);
        set_dependencies_gradient_term_RHS(var_index, string_gradn);
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

std::vector<scalarvalueType> eta_values(num_ops);
std::vector<scalarvalueType> dndt_values(num_ops);
std::vector<scalarvalueType> mu_values(num_muFields);
std::vector<scalargradType> mu_gradients(num_muFields);

for (unsigned int i=0; i<num_ops; ++i){
	eta_values[i] = variable_list.get_scalar_value(i);
	dndt_values[i] = variable_list.get_scalar_value(i+num_ops+num_muFields);
}
for (unsigned int i=0; i<num_muFields; ++i){
	mu_values[i] = variable_list.get_scalar_value(i+num_ops);
	mu_gradients[i] = variable_list.get_scalar_gradient(i+num_ops);
}
scalarvalueType sum_nsq = 0.0;
std::vector<scalarvalueType> sum_ops_sq_phase(num_phases, constV(0.0));
for (unsigned int i=0; i<num_ops; ++i){
	sum_nsq += eta_values[i]*eta_values[i];
    sum_ops_sq_phase[phase_index[i]] += eta_values[i]*eta_values[i];
}
std::vector<scalarvalueType> h(num_phases, constV(0.0));
std::vector<std::vector<scalarvalueType>>
    dhdn(num_phases, std::vector<scalarvalueType> (num_ops));

for (unsigned int phase_a=0; phase_a<num_phases; ++phase_a){
    for (unsigned int op=0; op<num_ops; ++op){
        scalarvalueType eta_sq = eta_values[op]*eta_values[op];
        if(phase_index[op] == phase_a){
            h[phase_a] += eta_sq/sum_nsq;
            dhdn[phase_a][op] = 2.0*eta_values[op]*(sum_nsq-sum_ops_sq_phase[phase_a])/(sum_nsq*sum_nsq);
        }
        else{
            dhdn[phase_a][op] = 2.0*eta_values[op]*(constV(0.0)-sum_ops_sq_phase[phase_a])/(sum_nsq*sum_nsq);
        }
    }
}

std::vector<scalarvalueType> dmudtValue(num_muFields);
std::vector<scalargradType> dmudtGrad(num_muFields);

for (unsigned int i=0; i<num_muFields; ++i){
    scalarvalueType susceptibility = 0.0;
    for (unsigned int j=0; j<num_phases; ++j){
        susceptibility += h[j]/(Va*Va*kWell[j][i]);
    }
    dmudtGrad[i] = -M*mu_gradients[i]/susceptibility;
    dmudtValue[i] = 0.0;
    for (unsigned int k=0; k<num_phases; ++k){
        for (unsigned int j=0; j<num_ops; ++j){
            scalarvalueType drhodn_part = dhdn[k][j]*(mu_values[i]/(Va*kWell[k][i]) + constV(cmin[k][i]))/Va;
            dmudtValue[i] -= dndt_values[j]*drhodn_part/susceptibility;
        }
    }
}

for (unsigned int i=0; i<num_ops; ++i){
    variable_list.set_scalar_value_term_RHS(i,eta_values[i] +
        dndt_values[i]*userInputs.dtValue);
}
for (unsigned int i=0; i<num_muFields; ++i){
    variable_list.set_scalar_value_term_RHS(i+num_ops,mu_values[i] +
        dmudtValue[i]*userInputs.dtValue);
    variable_list.set_scalar_gradient_term_RHS(i+num_ops,dmudtGrad[i]*userInputs.dtValue);
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

std::vector<scalarvalueType> eta_values(num_ops);
std::vector<scalargradType> eta_gradients(num_ops);
std::vector<scalarvalueType> mu_values(num_muFields);

for (unsigned int i=0; i<num_ops; ++i){
 	eta_values[i] = variable_list.get_scalar_value(i);
	eta_gradients[i] = variable_list.get_scalar_gradient(i);
}

for (unsigned int i=0; i<num_muFields; ++i){
	mu_values[i] = variable_list.get_scalar_value(i+num_ops);
}

std::vector<scalarvalueType> omegaC(num_phases);
scalarvalueType sum_nsq = 0.0;

for (unsigned int i=0; i<num_ops; ++i){
    sum_nsq += eta_values[i]*eta_values[i];
}

for (unsigned int i=0; i<num_phases; ++i){
    omegaC[i] = fWell[i];
    for (unsigned int j=0; j<num_muFields; ++j){
        omegaC[i] += -0.5*mu_values[j]*mu_values[j]/constV(Va*Va*kWell[i][j])
            - mu_values[j]*cmin[i][j]/Va;
    }
}

std::vector<scalarvalueType> dndtValue(num_ops);
std::vector<scalargradType> dndtGrad(num_ops);

for (unsigned int i=0; i < num_ops; ++i){
    dndtValue[i] = m0*(eta_values[i]*eta_values[i]*eta_values[i] - eta_values[i]);
    dndtGrad[i] = kappa*eta_gradients[i];
    dndtValue[i] += 2.0*eta_values[i]/sum_nsq * omegaC[phase_index[i]];
    for (unsigned int j=0; j<num_ops; ++j){//fix?
        if(i != j){
            dndtValue[i] += m0*2.0*eta_values[i]*gamma*eta_values[j]*eta_values[j];
        }
        dndtValue[i] -= 2.0*eta_values[j]*eta_values[j]*eta_values[i]/(sum_nsq*sum_nsq)
                        *omegaC[phase_index[j]];
    }
}

for (unsigned int i=0; i < num_ops; ++i){
    variable_list.set_scalar_value_term_RHS(i+num_ops+num_muFields,-L*dndtValue[i]);
    variable_list.set_scalar_gradient_term_RHS(i+num_ops+num_muFields,-L*dndtGrad[i]);
}

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
