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

    std::vector<double> c0 = userInputs.get_model_constant_double_array("c0");
    double center[3] = {0.5*userInputs.domain_size[0],0.5*userInputs.domain_size[1],0.0*userInputs.domain_size[2]};
    double pi = 3.141592653589793238;
    double x = p[0] - center[0];
    double y = p[1] - center[1];
    double z = p[2] - center[2];
    if(dim<3){z=0;}
    double intf = std::sqrt(0.5*m0/kappa);

    double r_A = 0.5*x + std::sqrt(3.0/4.0)*y;
    double r_B = 0.5*x - std::sqrt(3.0/4.0)*y;
    double r_C = -x;

    double power = userInputs.get_model_constant_double("power");
    double spacing = userInputs.get_model_constant_double("spacing");
    double scaled_r = std::pow(std::pow(r_A*r_A,power) + std::pow(r_B*r_B,power) + std::pow(r_C*r_C,power), 0.5/power);
    double wave = std::cos((2.0*pi*scaled_r/spacing) + std::atan2(x,y))*spacing*0.5/pi;

    double tanh_profile_z = 0.5*(1 - std::tanh(intf*(z-0.5)));
    
    op_vals[0] = 1.0 - tanh_profile_z;//liquid
    op_vals[1] = tanh_profile_z * 0.5*(1.0 - std::tanh(intf*wave));
    op_vals[2] = tanh_profile_z * 0.5*(1.0 + std::tanh(intf*wave));
    //mu_vals[0] = 12.0;//component 0
    for(unsigned int mu_index = 0; mu_index<num_muFields; ++mu_index){
        mu_vals[0] = kWell[0][mu_index]*(c0[mu_index]-cmin[0][mu_index]);
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
