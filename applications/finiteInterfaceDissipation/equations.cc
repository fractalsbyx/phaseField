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
    const unsigned int num_comps{1};
    const unsigned int num_phases{1};
    std::string string_val_n = "";
    std::string string_val_c = "";
    std::string string_grad_n = "";
    std::string string_grad_c = "";
    for (unsigned int op_index=0; op_index<num_ops; op_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(op_index));
        string_val_n.append(var_name+",");
        string_grad_n.append("grad("+var_name+"),");
        set_variable_name				(op_index, var_name);
    	set_variable_type				(op_index, SCALAR);
    	set_variable_equation_type		(op_index, EXPLICIT_TIME_DEPENDENT);
        set_allowed_to_nucleate			(op_index, (op_index>0));
        set_need_value_nucleation		(op_index, true);
    }
    for (unsigned int phase_index=0; phase_index<num_phases; phase_index++){
        for (unsigned int comp_index=0; comp_index<num_comps; comp_index++){
            std::string var_name = "c_";
            var_name.append(std::to_string(phase_index + "_"));
            var_name.append(std::to_string(comp_index));
            string_val_c.append(var_name+",");
            string_grad_c.append("grad("+var_name+"),");
            set_variable_name				(num_ops+(num_comps*phase_index)+comp_index, var_name);
            set_variable_type				(num_ops+(num_comps*phase_index)+comp_index, SCALAR);
            set_variable_equation_type		(num_ops+(num_comps*phase_index)+comp_index, EXPLICIT_TIME_DEPENDENT);
            set_need_value_nucleation		(num_ops+(num_comps*phase_index)+comp_index, true);
        }
    }
    std::cout << string_val_n << " | " << string_val_c << " | "
            << string_grad_n << " | " << string_grad_c << "\n";
    std::string dep_val_n = string_val_n+string_val_c;
    std::string dep_val_c = string_val_n+string_val_c;
    std::string dep_grad_n = string_val_n+string_grad_n;
    std::string dep_grad_c = string_val_n+string_val_c+string_grad_c;
    dep_val_n.pop_back();
    dep_val_c.pop_back();
    for (unsigned int op_index=0; op_index<num_ops; op_index++){
        set_dependencies_value_term_RHS(op_index, dep_val_n);
        set_dependencies_gradient_term_RHS(op_index, dep_grad_n);
    }
    for (unsigned int phase_index = 0; phase_index < num_phases; phase_index++){
        for (unsigned int comp_index = 0; comp_index < num_comps; comp_index++){
            set_dependencies_value_term_RHS(num_ops+(num_comps*phase_index)+comp_index, dep_val_c);
            set_dependencies_gradient_term_RHS(num_ops+(num_comps*phase_index)+comp_index, dep_grad_c);
        }
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
// "Order Parameter" for each phase [phase_name]
std::unordered_map<std::string, scalarvalueType> phi_phase_val(num_phases);
std::unordered_map<std::string, scalargradType> phi_phase_grad(num_phases);
// Local composition [component_name]
std::unordered_map<std::string, scalarvalueType> c_val;
std::unordered_map<std::string, scalargradType> c_grad;
// Composition within phase [phase_name][component_name]
std::unordered_map<std::string, std::unordered_map<std::string, scalarvalueType>> c_phase_val;
std::unordered_map<std::string, std::unordered_map<std::string, scalargradType>> c_phase_grad;
// Sublattice occupation [phase_name][sublattice index s][component name]
std::unordered_map<std::string, std::vector<std::unordered_map<std::string, scalarvalueType>>> y_val;
std::unordered_map<std::string, std::vector<std::unordered_map<std::string, scalargradType>>> y_grad;
//
std::unordered_map<std::string, std::vector<std::unordered_map<std::string, scalarvalueType>>> dydt_val;
std::unordered_map<std::string, std::vector<std::unordered_map<std::string, scalargradType>>> dydt_grad;
// something
std::unordered_map<std::string, std::vector<scalarvalueType>> phi_dcdt_val(num_phases, std::vector<scalarvalueType>(num_comps));
std::unordered_map<std::string, std::vector<scalargradType>> phi_dcdt_grad(num_phases, std::vector<scalargradType>(num_comps));

std::vector<scalarvalueType> deltaG(num_phases, constV(0.0));

// Retrieve fields
unsigned int var_index = 0;
for (auto phase = phase_names.begin(); phase < phase_names.end(); phase++){
    for (unsigned int s=0; s<num_sublattices[*phase]; s++){
        if (sublattice_comps[*phase][s].size()>1){ // saving memory for fixed sublattices
        for (auto comp = sublattice_comps[*phase][s].begin(); comp<sublattice_comps[*phase][s].end(); comp++){
            y_val[*phase][s][*comp] = variable_list.get_scalar_value(var_index);
            y_grad[*phase][s][*comp] = variable_list.get_scalar_gradient(var_index++);
            c_phase_val[*phase][*comp] += num_sites[*phase][s] * y_val[*phase][s][*comp] / total_num_sites[*phase];
            c_phase_grad[*phase][*comp] += num_sites[*phase][s] * y_grad[*phase][s][*comp] / total_num_sites[*phase];
        }
        }
        else{
            auto comp = sublattice_comps[*phase][s].begin();
            c_phase_val[*phase][*comp] += num_sites[*phase][s] * 1.0 / total_num_sites[*phase];
            c_phase_grad[*phase][*comp] += num_sites[*phase][s] * 1.0 / total_num_sites[*phase];
        }
    }
}
for (unsigned int i=0; i<num_ops; ++i){
	phi_val[i] =  variable_list.get_scalar_value(var_index);
	phi_grad[i] = variable_list.get_scalar_gradient(var_index++);
    auto phase = op_phase_name.begin()+i; // Could just use index, want to be consistent with rest of code.
    phi_phase_val[*phase] += phi_val[i];   // This and phase comps need order parameters to sum to 1,
    phi_phase_grad[*phase] += phi_grad[i]; // otherwise this must be done later with interpolation function.
    for (auto comp = comp_names.begin(); comp<comp_names.end(); comp++){
        c_val[*comp] += c_phase_val[*phase][*comp] * phi_val[i];
        c_grad[*comp] +=  c_phase_grad[*phase][*comp] * phi_val[i] + c_phase_val[*phase][*comp] * phi_grad[i];
    }
}

//
for (unsigned int j=0; j<num_phases; j++){
    for (unsigned int k=0; k<num_comps; k++){
        phi_dcdt_grad[j][k] = phi_phase[j]*D[j]*c_grad[j][k];
        phi_dcdt_val[j][k] = phi_phase[j]
            *(P*(1-phi_phase[j])*(sum_mu[k] - mu[j][k])
            -phi_phase[j]*sum_dphidt_c[j]);
    }
}

for (unsigned int i=0; i<num_ops; ++i){
    dphidt_val[i] = sigma*(pi/l_gb)*(pi/l_gb)*(phi_val[i]-0.5);
    dphidt_grad[i] = sigma*phi_grad[i];
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

std::vector<scalarvalueType> phi_values(num_ops);
std::vector<scalargradType> phi_gradients(num_ops);
std::vector<scalarvalueType> c_values(num_cFields);

for (unsigned int i=0; i<num_ops; ++i){
 	phi_values[i] = variable_list.get_scalar_value(i);
	phi_gradients[i] = variable_list.get_scalar_gradient(i);
}

for (unsigned int i=0; i<num_cFields; ++i){
	c_values[i] = variable_list.get_scalar_value(i+num_ops);
}

std::vector<scalarvalueType> omegaC(num_phases);
scalarvalueType sum_nsq = 0.0;

for (unsigned int i=0; i<num_ops; ++i){
    sum_nsq += phi_values[i]*phi_values[i];
}

for (unsigned int i=0; i<num_phases; ++i){
    omegaC[i] = fWell[i];
    for (unsigned int j=0; j<num_cFields; ++j){
        omegaC[i] += -0.5*c_values[j]*c_values[j]/constV(Va*Va*kWell[i][j])
            - c_values[j]*cmin[i][j]/Va;
    }
}

std::vector<scalarvalueType> dndtValue(num_ops);
std::vector<scalargradType> dndtGrad(num_ops);

for (unsigned int i=0; i < num_ops; ++i){
    dndtValue[i] = m0*(phi_values[i]*phi_values[i]*phi_values[i] - phi_values[i]);
    dndtGrad[i] = kappa*phi_gradients[i];
    dndtValue[i] += 2.0*phi_values[i]/sum_nsq * omegaC[phase_index[i]];
    for (unsigned int j=0; j<num_ops; ++j){//fix?
        if(i != j){
            dndtValue[i] += m0*2.0*phi_values[i]*gamma*phi_values[j]*phi_values[j];
        }
        dndtValue[i] -= 2.0*phi_values[j]*phi_values[j]*phi_values[i]/(sum_nsq*sum_nsq)
                        *omegaC[phase_index[j]];
    }
}

//std::vector<double> forcingTerm(num_ops, 0.0);
std::vector<dealii::VectorizedArray<double>> forcingTerm(num_ops, constV(0.0));
std::vector<dealii::VectorizedArray<double>> mob_term(num_ops, constV(1.0));
double interface_width = 1.0/std::sqrt(0.5*m0/kappa);
seedNucleus(q_point_loc, forcingTerm, mob_term, interface_width);

for (unsigned int i=0; i < num_ops; ++i){
    //variable_list.set_scalar_value_term_RHS(i+num_ops+num_cFields,-L*(dndtValue[i])*mob_term[i]+forcingTerm[i]);
    //variable_list.set_scalar_gradient_term_RHS(i+num_ops+num_cFields,-L*(dndtGrad[i])*mob_term[i]);
    variable_list.set_scalar_value_term_RHS(i+num_ops+num_cFields,-L*(dndtValue[i]-nuc_force*forcingTerm[i]));
    variable_list.set_scalar_gradient_term_RHS(i+num_ops+num_cFields,-L*(dndtGrad[i]));
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

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	std::vector<dealii::VectorizedArray<double>> & source_term,
	std::vector<dealii::VectorizedArray<double>> & mob_term,
    double interface_coeff) const {

    for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){
        if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){
            // Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
            unsigned int nuc_op_index = thisNucleus->orderParameterIndex;
            dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(nuc_op_index), q_point_loc, nuc_op_index);

            for (unsigned i=0; i<mob_term[nuc_op_index].size();i++){
                if (weighted_dist[i] <= 1.0){
                    mob_term[nuc_op_index][i] = 0.0;

                    // Seed a nucleus if it was added to the list of nuclei this time step
                    if (true/*thisNucleus->seedingTimestep == this->currentIncrement*/){
                        // Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
                        dealii::Point<dim,double> q_point_loc_element;
                        for (unsigned int j=0; j<dim; j++){
                            q_point_loc_element(j) = q_point_loc(j)[i];
                        }
                        double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex), q_point_loc_element, thisNucleus->orderParameterIndex);

                        double avg_semiaxis = 0.0;
                        for (unsigned int j=0; j<dim; j++){
                            avg_semiaxis += thisNucleus->semiaxes[j];
                        }
                        avg_semiaxis /= dim;

                        source_term[nuc_op_index][i] += 0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));
                    }
                }
            }
        }
    }
}