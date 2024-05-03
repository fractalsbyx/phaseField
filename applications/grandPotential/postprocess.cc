// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){
    const unsigned int num_ops{2};
    const unsigned int num_muFields{2};
    std::string string_valn = "";
    std::string string_valmu = "";
for (unsigned int n_index=0; n_index<num_ops; n_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(n_index));
        string_valn.append(var_name+",");      
    }
    for (unsigned int c_index=0; c_index<num_muFields; c_index++){
        std::string mu_name = "mu";
		std::string c_name = "c";
        mu_name.append(std::to_string(c_index));
		c_name.append(std::to_string(c_index));
        string_valmu.append(mu_name+",");
        set_variable_name(c_index,c_name);
    	set_variable_type(c_index,SCALAR);
		set_output_integral(c_index, true);
    }
	std::string dep = string_valmu + string_valn;
	dep.pop_back();
	for (unsigned int c_index=0; c_index<num_muFields; c_index++){
		set_dependencies_value_term_RHS(c_index, dep);
	}

	set_variable_name(num_muFields+0, "sum_nsq");
    set_variable_type(num_muFields+0,SCALAR);
	set_dependencies_value_term_RHS(num_muFields+0, dep);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.h. It takes in "variable_list" and "q_point_loc" as inputs and outputs two terms in
// the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.h' for assigning the values/derivatives of
// the primary variables.

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    // --- Getting the values and derivatives of the model variables ---
std::vector<scalarvalueType> eta_values(num_ops);
std::vector<scalarvalueType> mu_values(num_muFields);
std::vector<scalarvalueType> c_values(num_muFields);

for (unsigned int i=0; i<num_ops; ++i){
 	eta_values[i] = variable_list.get_scalar_value(i);
}

for (unsigned int i=0; i<num_muFields; ++i){
	mu_values[i] = variable_list.get_scalar_value(i+num_ops);
}

scalarvalueType sum_nsq = 0.0;

for (unsigned int i=0; i<num_ops; ++i){
    sum_nsq += eta_values[i]*eta_values[i];
}

for (unsigned int op_index=0; op_index<num_ops; ++op_index){
	scalarvalueType h = eta_values[op_index]*eta_values[op_index]/sum_nsq;
	for (unsigned int mu_index=0; mu_index<num_muFields; ++mu_index){
    	c_values[mu_index] += h*(mu_values[mu_index]/(Va*kWell[phase_index[op_index]][mu_index]) + cmin[phase_index[op_index]][mu_index]);
	}
}

for (unsigned int i=0; i<num_muFields; ++i){
    pp_variable_list.set_scalar_value_term_RHS(i, c_values[i]);
}
pp_variable_list.set_scalar_value_term_RHS(num_muFields+0, sum_nsq);
}
