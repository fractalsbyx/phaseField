// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index
    unsigned int N = userInputs.number_of_variables-1;
    scalar_IC = 0;



    //if (index > 37 && index < 76)
        //scalar_IC = rho_i[index-38] + 2.0*icamplitude*(dist(rng)-0.5);


    if (index == N)
      scalar_IC = 0;


    /*if (index == 38)
        scalar_IC = rho_i[0];
    if (index == 39)
        scalar_IC = rho_i[1];
    if (index == 40)
        scalar_IC = rho_i[2];*/
    //if (index == N+1)
        //scalar_IC = 0;
	// --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
