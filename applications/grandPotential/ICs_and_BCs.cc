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

	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".
    /*
    double center[1][3] = {{0.5,0.5,0.5}};
    double rad[1] = {2.0};
    double distance;
    scalar_IC = 0.0;
    // Initial condition for the concentration field
    if ((index == 0) or (index == 1)){
        for (unsigned int i=0; i<1; i++){
            distance = 0.0;
            for (unsigned int dir = 0; dir < dim; dir++){
                distance += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*
                    (p[dir]-center[i][dir]*userInputs.domain_size[dir]);
            }
            distance = std::sqrt(distance);
            if(index == 0){
            //    scalar_IC = 1.0;
                scalar_IC = 0.5*(1+std::tanh(2.0*(distance-rad[i])));
            }
            if(index == 1){
                scalar_IC = 1.0 - 0.5*(1.0+std::tanh(2.0*(distance-rad[i])));
            }
        }
    }
    if (index > 3){
        //scalar_IC = 0.1*(dist(rng)-0.05);
        scalar_IC = 0;
    }*/
    mu_ini = std::vector<std::vector<double>> (num_phases, std::vector<double> (num_comps));
    for(unsigned int phase = 0; phase < num_phases; phase++){
        for(unsigned int comp = 0; comp < num_comps; comp++){
            mu_ini[phase][comp] = 0.0;//kWell[phase][comp]*(c0[comp] - cmin[phase][comp]);
        }
    }

    double center[3] = {5.0,5.0,5.0};
    double rad = 2.0;
    double dist2 = 0.0;
    double dist = 0.0;
    for(unsigned int xyz = 0; xyz<dim; xyz++){
        dist2 += (p[xyz]-center[xyz])*(p[xyz]-center[xyz]);
    }
    dist = std::sqrt(dist2);
    double tanh_profile = 0.5*(1+std::tanh(2.0*(dist-rad)));
    if (index==0){
        scalar_IC = tanh_profile;
    }
    else if (index==1){
        scalar_IC = 1 - tanh_profile;
    }
    else if (index<=num_phases){
        scalar_IC = 0;
    }
    else if (index<=num_phases+num_comps){
        scalar_IC = tanh_profile*mu_ini[0][index-num_phases] + (1-tanh_profile)*mu_ini[1][index-num_phases];
    }
    else{
        scalar_IC = 0.0;
    }
    // Initial condition for the order parameter field

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
