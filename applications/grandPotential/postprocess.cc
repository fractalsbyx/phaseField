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
	// Variable 0
    set_variable_name                (0,"sum_nsq");
    set_variable_type                (0,SCALAR);

    set_dependencies_value_term_RHS(0, "n0, n1, n2");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral             (0,false);

    // Variable 1
    set_variable_name                (1,"gb");
    set_variable_type                (1,SCALAR);

    set_dependencies_value_term_RHS(1, "n0, n1, n2");
    set_dependencies_gradient_term_RHS(1, "");

    set_output_integral             (1,false);
	// Var 2
	set_variable_name                (2,"sum_n");
    set_variable_type                (2,SCALAR);

    set_dependencies_value_term_RHS(2, "n0, n1, n2");
    set_dependencies_gradient_term_RHS(2, "");

    set_output_integral             (2,false);
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
    scalarvalueType sum_nsq = constV(0.0);
	scalarvalueType sum_n   = constV(0.0);
    for (unsigned int i=0; i<num_phases; i++){
        scalarvalueType ni = variable_list.get_scalar_value(i);
        sum_nsq += ni*ni;
		sum_n += ni;
    }
	scalarvalueType gb = constV(1.0) - sum_nsq;
	// The concentration and its derivatives

	// Submit
	pp_variable_list.set_scalar_value_term_RHS(0, sum_nsq);
	pp_variable_list.set_scalar_value_term_RHS(1, gb);
	pp_variable_list.set_scalar_value_term_RHS(2, sum_n);
}
