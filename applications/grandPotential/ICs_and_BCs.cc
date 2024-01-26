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

    // std::vector<std::vector<double>> mu_ini(num_phases, std::vector<double> (num_comps));
    // for(unsigned int phase = 0; phase < num_phases; phase++){
    //     for(unsigned int comp = 0; comp < num_muFields; comp++){
    //         mu_ini[phase][comp] = 0.0;//kWell[phase][comp]*(c0[comp] - cmin[phase][comp]);
    //     }
    // }



    double center[3] = {5.0,5.0,5.0};
    double rad = 2.0;
    double xwidth = 2.5;
    double dist2 = 0.0;
    double dist = 0.0;
    double xdist = p[0]-center[0];
    for(unsigned int xyz = 0; xyz<dim; xyz++){
        dist2 += (p[xyz]-center[xyz])*(p[xyz]-center[xyz]);
    }
    dist = std::sqrt(dist2);
    double tanh_profile = 0.5*(1+std::tanh(-2.0*(dist-rad)));
    double tanh_profileB= 0.5*(1+std::tanh(-2.0*(std::abs(xdist)-xwidth)));
    if (index<=num_ops){
        if (index==0){
            scalar_IC = (1 - tanh_profile)*(1-tanh_profileB);
        }
        else if (index==1){
            scalar_IC = tanh_profile;
        }
        else if (index==2){
            scalar_IC = tanh_profileB*(1-tanh_profile);
        }
        else{
            scalar_IC = 0.0;
        }
    }
    else if (index<=num_ops+num_comps){
        scalar_IC = 0;
        // for(unsigned int phase = 0; phase<num_phases; ++phase){
        //     scalar_IC += 0.0;// tanh_profile*mu_ini[0][index-num_phases] + (1-tanh_profile)*mu_ini[1][index-num_phases];
        // }
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
