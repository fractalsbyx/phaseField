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

    // Precalculating everything makes writing initial conditions easier. May take slightly more runtime.
    std::vector<double> op_vals(num_ops, 0.0);
    std::vector<double> mu_vals(num_muFields, 0.0);
    std::vector<double> h(num_ops, 0.0);

    // Custom coordinate system
    double center[3] = {0.5*userInputs.domain_size[0],0.5*userInputs.domain_size[1],0.0*userInputs.domain_size[2]};
    double x = p[0] - center[0];
    double y = p[1] - center[1];
    double z = p[2] - center[2];
    if(dim<3){z=0;}
    double r2 = x*x+y*y+z*z;

    // constant definitions
    double pi = 3.141592653589793238;
    double intf = std::sqrt(0.5*m0/kappa);
    std::vector<std::vector<double>> c0  = reshapeVector(userInputs.get_model_constant_double_array("c0"),
                                                                  num_phases, num_muFields);
    double rad = userInputs.get_model_constant_double("r0");

    // profile for making order parameters
    double tanh_circle = 0.5*(1.0 - std::tanh(0.5*intf*(r2-rad*rad)/rad));
    
    // make order parameters
    op_vals[0] = 1.0 - tanh_circle;//liquid
    op_vals[1] = tanh_circle;

    // Interpolation fields
    double sum_nsq = 0;
    for(unsigned int op_index = 0; op_index<num_ops; ++op_index){
        sum_nsq += op_vals[op_index]*op_vals[op_index];
    }
    for(unsigned int op_index = 0; op_index<num_ops; ++op_index){
        h[op_index] = op_vals[op_index]*op_vals[op_index]/sum_nsq;
    }

    // make mu fields
    for(unsigned int mu_index = 0; mu_index<num_muFields; ++mu_index){
        for(unsigned int op_index = 0; op_index<num_ops; ++op_index){
            mu_vals[mu_index] += h[op_index]
                                *Va*kWell[phase_index[op_index]][mu_index]
                                *(c0[phase_index[op_index]][mu_index]
                                    -cmin[phase_index[op_index]][mu_index]);
        }
        //mu_vals[mu_index] = Va*kWell[phase_index[0]][mu_index]*(c0[mu_index]-cmin[phase_index[0]][mu_index]);
    }

    // ===========================================================================
    // Submit fields
    // ===========================================================================
    scalar_IC = 0.0;
    for(unsigned int op_index = 0; op_index<num_ops; ++op_index){
        if(index==op_index){scalar_IC = op_vals[op_index];}
    }
    for(unsigned int mu_index = 0; mu_index<num_muFields; ++mu_index){
        if(index==num_ops+mu_index){scalar_IC = mu_vals[mu_index];}
    }
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
