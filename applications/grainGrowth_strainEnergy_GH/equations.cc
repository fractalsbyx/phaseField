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

    // GRAIN ORDER PARAMETERS
    unsigned int N=25;
    for (unsigned int var_index=0; var_index<N; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));

        std::string grad_var_name = "grad(n";
        grad_var_name.append(std::to_string(var_index));
        grad_var_name.append(")");

        set_variable_name                (var_index,var_name);
        set_variable_type                (var_index,SCALAR);
        set_variable_equation_type       (var_index,EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(var_index, var_name);
        set_dependencies_gradient_term_RHS(var_index, grad_var_name);
    }

    // DISLOCATION ORDER PARAMETERS
    unsigned int N_rho=25;
    for (unsigned int index=0; index<N_rho; index++){
        std::string var_nam = "rho";
        var_nam.append(std::to_string(index));


        set_variable_name                (index+N,var_nam);
        set_variable_type                (index+N,SCALAR);
        set_variable_equation_type       (index+N,EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(index+N, var_nam);
    }

    set_variable_name                (N+N_rho,"rho");
    set_variable_type                (N+N_rho,SCALAR);
    set_variable_equation_type       (N+N_rho,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(N+N_rho, "rho");

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

    unsigned int N = 25;
    unsigned int N_rho = 25;

    scalarvalueType fnV = constV(0.0);
    scalarvalueType fsV = constV(0.0);
    scalarvalueType sum_nsq = constV(0.0);
    scalarvalueType rhocalc = constV(0.0);
    scalarvalueType ni, nj, rhoi;//, rhocalc;
    scalargradType nix;
    
    std::vector<scalarvalueType> value_terms;
    value_terms.resize(N);
    std::vector<scalargradType> gradient_terms;
    gradient_terms.resize(N);

    //float t_n;
    //t_n = this->currentTime;

    // ------------------------------------------------------------
    // Perform calculations
    // ------------------------------------------------------------

 
    // CALCULATE SUM OF ORDER PARAMETERS SQUARED
    for (unsigned int i=0; i<N; i++){
        ni = variable_list.get_scalar_value(i);
        sum_nsq += ni*ni;
    }
    // INTERPOLATION FOR RHO FIELD
    for (unsigned int i=0; i<N; i++){
        ni = variable_list.get_scalar_value(i);
        rhoi = variable_list.get_scalar_value(i+N);
        rhocalc += ni*ni*rhoi/sum_nsq;
    }
    // CALCULATE FREE ENERGY TERMS
    for (unsigned int i=0; i<N; i++){
        ni = variable_list.get_scalar_value(i);
        nix = variable_list.get_scalar_gradient(i);
        rhoi = variable_list.get_scalar_value(i+N);
        fnV = - ni + ni*ni*ni;
        for (unsigned int j=0; j<N; j++){
            if (i != j){
                nj = variable_list.get_scalar_value(j);
                fnV += constV(2.0*alpha) * ni * nj*nj;
            }
        }
        fnV = m0*fnV;
        //fsV = 0.0; // uncomment to turn off stored-energy
        fsV = with_stored*constV(2.0*fs_const)*ni*(rhoi-rhocalc)/sum_nsq; // Guanglong's expression for dp/dn
        value_terms[i] = ni-constV(userInputs.dtValue*MnV)*(fnV + fsV);
        gradient_terms[i] = constV(-userInputs.dtValue*KnV*MnV)*nix;
    }
    // SUBMIT TERMS FOR GOVERNING EQUATIONS
    for (unsigned int i=0; i<N; i++){
        variable_list.set_scalar_value_term_RHS(i,value_terms[i]);
        variable_list.set_scalar_gradient_term_RHS(i,gradient_terms[i]);
    }
    for (unsigned int i=0; i < N_rho; i++){
      rhoi = variable_list.get_scalar_value(i+N);
      variable_list.set_scalar_value_term_RHS(i+N,rhoi);
      variable_list.set_scalar_gradient_term_RHS(i+N,nix*constV(0.0));
    }

    variable_list.set_scalar_value_term_RHS(N+N_rho,rhocalc);

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

//t_n = this->currentTime;
    //if (t_n < 0.01)
    /*for (unsigned int i=0; i < N; i++){
      rhoi = variable_list.get_scalar_value(i+N);
      ni = variable_list.get_scalar_value(i);
      variable_list.set_scalar_value_term_RHS(i+N,ni);//*rho_i[i]);
    }*/
