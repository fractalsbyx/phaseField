// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain.
#include <random>


void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"feature_ids");
	set_variable_type				(0,SCALAR);

    // For the input file 'parameters.in'
    set_dependencies_value_term_RHS(0, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19");

    // For the input file 'parameters_large_2D.in'
    //set_dependencies_value_term_RHS(0, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11");

    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,false);

    // Variable 1
    set_variable_name				(1,"op_ids");
	set_variable_type				(1,SCALAR);

    // For the input file 'parameters.in'
    set_dependencies_value_term_RHS(1, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19");

    // For the input file 'parameters_large_2D.in'
    //set_dependencies_value_term_RHS(1, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11");

    set_dependencies_gradient_term_RHS(1, "");

    set_output_integral         	(1,false);


    // Dislocation Densities
    set_variable_name               (2,"rho0");
    set_variable_type               (2,SCALAR);

    set_dependencies_value_term_RHS(2, "n0");
    set_dependencies_gradient_term_RHS(2, "");

    set_output_integral             (2,false);

    set_variable_name               (3,"rho1");
    set_variable_type               (3,SCALAR);

    set_dependencies_value_term_RHS(3, "n1");
    set_dependencies_gradient_term_RHS(3, "");

    set_output_integral             (3,false);

    set_variable_name               (4,"rho2");
    set_variable_type               (4,SCALAR);

    set_dependencies_value_term_RHS(4, "n2");
    set_dependencies_gradient_term_RHS(4, "");

    set_output_integral             (4,false);

    set_variable_name               (5,"rho3");
    set_variable_type               (5,SCALAR);

    set_dependencies_value_term_RHS(5, "n3");
    set_dependencies_gradient_term_RHS(5, "");

    set_output_integral             (5,false);

    set_variable_name               (6,"rho4");
    set_variable_type               (6,SCALAR);

    set_dependencies_value_term_RHS(6, "n4");
    set_dependencies_gradient_term_RHS(6, "");

    set_output_integral             (6,false);

    set_variable_name               (7,"rho5");
    set_variable_type               (7,SCALAR);

    set_dependencies_value_term_RHS(7, "n5");
    set_dependencies_gradient_term_RHS(7, "");

    set_output_integral             (7,false);

    set_variable_name               (8,"rho6");
    set_variable_type               (8,SCALAR);

    set_dependencies_value_term_RHS(8, "n6");
    set_dependencies_gradient_term_RHS(8, "");

    set_output_integral             (8,false);

    set_variable_name               (9,"rho7");
    set_variable_type               (9,SCALAR);

    set_dependencies_value_term_RHS(9, "n7");
    set_dependencies_gradient_term_RHS(9, "");

    set_output_integral             (9,false);

    set_variable_name               (10,"rho8");
    set_variable_type               (10,SCALAR);

    set_dependencies_value_term_RHS(10, "n8");
    set_dependencies_gradient_term_RHS(10, "");

    set_output_integral             (10,false);

    set_variable_name               (11,"rho9");
    set_variable_type               (11,SCALAR);

    set_dependencies_value_term_RHS(11, "n9");
    set_dependencies_gradient_term_RHS(11, "");

    set_output_integral             (11,false);

    set_variable_name               (12,"rho10");
    set_variable_type               (12,SCALAR);

    set_dependencies_value_term_RHS(12, "n10");
    set_dependencies_gradient_term_RHS(12, "");

    set_output_integral             (12,false);

    set_variable_name               (13,"rho11");
    set_variable_type               (13,SCALAR);

    set_dependencies_value_term_RHS(13, "n11");
    set_dependencies_gradient_term_RHS(13, "");

    set_output_integral             (13,false);

    set_variable_name               (14,"rho12");
    set_variable_type               (14,SCALAR);

    set_dependencies_value_term_RHS(14, "n12");
    set_dependencies_gradient_term_RHS(14, "");

    set_output_integral             (14,false);

    set_variable_name               (15,"rho13");
    set_variable_type               (15,SCALAR);

    set_dependencies_value_term_RHS(15, "n13");
    set_dependencies_gradient_term_RHS(15, "");

    set_output_integral             (15,false);

    set_variable_name               (16,"rho14");
    set_variable_type               (16,SCALAR);

    set_dependencies_value_term_RHS(16, "n14");
    set_dependencies_gradient_term_RHS(16, "");

    set_output_integral             (16,false);

    set_variable_name               (17,"rho15");
    set_variable_type               (17,SCALAR);

    set_dependencies_value_term_RHS(17, "n15");
    set_dependencies_gradient_term_RHS(17, "");

    set_output_integral             (17,false);

    set_variable_name               (18,"rho16");
    set_variable_type               (18,SCALAR);

    set_dependencies_value_term_RHS(18, "n16");
    set_dependencies_gradient_term_RHS(18, "");

    set_output_integral             (18,false);

    set_variable_name               (19,"rho17");
    set_variable_type               (19,SCALAR);

    set_dependencies_value_term_RHS(19, "n17");
    set_dependencies_gradient_term_RHS(19, "");

    set_output_integral             (19,false);

    set_variable_name               (20,"rho18");
    set_variable_type               (20,SCALAR);

    set_dependencies_value_term_RHS(20, "n18");
    set_dependencies_gradient_term_RHS(20, "");

    set_output_integral             (20,false);

    set_variable_name               (21,"rho19");
    set_variable_type               (21,SCALAR);

    set_dependencies_value_term_RHS(21, "n19");
    set_dependencies_gradient_term_RHS(21, "");

    set_output_integral             (21,false);

    set_variable_name               (22,"rhocalc");
    set_variable_type               (22,SCALAR);

    set_dependencies_value_term_RHS(22, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19");
    set_dependencies_gradient_term_RHS(22, "");

    set_output_integral             (22,false);

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

scalarvalueType ni, n0;

scalarvalueType max_val = constV(-100.0);
scalarvalueType max_op = constV(100.0);
scalarvalueType Rho = constV(0.0);
scalarvalueType sum_nsq = constV(0.0);
scalarvalueType rhocalc = constV(0.0);

for (unsigned int i=0; i<userInputs.number_of_variables; i++){
    ni = variable_list.get_scalar_value(i);

    for (unsigned int v=0; v<ni.size();v++){
        if (ni[v] > max_val[v]){
            max_val[v] = ni[v];
            max_op[v] = i;
        }
    }
}

// Calculate sum of eta squared for Rho interpolation
for (unsigned int i=0; i<userInputs.number_of_variables; i++){
    ni = variable_list.get_scalar_value(i);
    sum_nsq += ni*ni;
}


// Step 1: Get list of grain centers (for 20 variables)
std::vector<std::vector<dealii::Point<dim>>> grain_center_list(20);
std::vector<std::vector<unsigned int>> grain_center_ID(20);

// Loop over grains
for (unsigned int g=0; g<this->simplified_grain_representations.size(); g++){
    
    unsigned int OP = this->simplified_grain_representations[g].getOrderParameterId();

    // Loop over order parameters and push back center and ID to list
    for (unsigned int i=0; i<userInputs.number_of_variables; i++){
        if (OP == i){
            // Step 1: Generate the center position coordinates for each grain
            dealii::Point<dim> grain_center;
            for (unsigned int d=0;d<dim;d++){
                grain_center(d) = this->simplified_grain_representations[g].getCenter()(d);
            }
            grain_center_list[i].push_back(grain_center);
            grain_center_ID[i].push_back(this->simplified_grain_representations[g].getGrainId());
        }//endif OP
    }
}

// In this section, change the values of the array grain_center_ID from FeatureIds to dislocation density
// The Feature IDs are reordered to 0, 1,..., 121 where 0 and 1 are unindexed and outer regions. Feature IDs
// 2 through 121 correspond to the MATLAB data array GNDlistPink
std::vector<std::vector<float>> grain_center_rho(20);



// Step 2: Loop over domain and calculate distance to grain centers at each point
// Assign the value of the grain to the current point.

// Loop over each order parameter (0 to 19) and calculate the 
for (unsigned int i=0; i<userInputs.number_of_variables; i++){

    for (unsigned int v=0; v<ni.size();v++){

        // Get the position vector of the current element
        dealii::Point<dim> q_point_loc_nonvec;
        for (unsigned int d=0;d<dim;d++){
            q_point_loc_nonvec(d) = q_point_loc(d)[v];
        }

        double closest_grain;
        double closest_grain_dist = 500.0;
        // Loop over grain centers and find the smallest distance
        for (unsigned int grain=0; grain < grain_center_list[i].size(); grain++){

            // Calculate the distance from the center of each grain
            double dist = 0.0;
            for (unsigned int d=0;d<dim;d++){
                dist += (q_point_loc_nonvec(d) - grain_center_list[i][grain](d))*(q_point_loc_nonvec(d) - grain_center_list[i][grain](d));
            }
            dist = std::sqrt(dist);

            // If the distance to the current grain is smaller, set it as the closet grain
            if (dist < closest_grain_dist){
                closest_grain = grain_center_ID[i][grain];
                closest_grain_dist = dist;
            }
        }//endfor grain

        // Set the value of Rho equal to the ID of the closest grain
        Rho[v] = double(closest_grain);

    }//endfor v

    ni = variable_list.get_scalar_value(i);
    ni = std::max(ni, constV(0.0));

    Rho = Rho*ni;
    Rho = std::max(Rho, constV(0.0));

    // Interpolation for rho field
    rhocalc += ni*ni*Rho/sum_nsq;

    pp_variable_list.set_scalar_value_term_RHS(i+2, Rho);
}

pp_variable_list.set_scalar_value_term_RHS(22, rhocalc);

scalarvalueType feature_ids = constV(-1.0);
for (unsigned int v=0; v<ni.size();v++){

    for (unsigned int g=0; g<this->simplified_grain_representations.size(); g++){

        unsigned int max_op_nonvec = (unsigned int)std::abs(max_op[v]);
        if (this->simplified_grain_representations[g].getOrderParameterId() == max_op_nonvec){

            // Get the position vector of the current element
            dealii::Point<dim> q_point_loc_nonvec;
            for (unsigned int d=0;d<dim;d++){
                q_point_loc_nonvec(d) = q_point_loc(d)[v];
            }

            // Calculate the distance from the center of grain g
            double dist = 0.0;
            for (unsigned int d=0;d<dim;d++){
                dist += (q_point_loc_nonvec(d)-this->simplified_grain_representations[g].getCenter()(d))*(q_point_loc_nonvec(d)-this->simplified_grain_representations[g].getCenter()(d));
            }
            dist = std::sqrt(dist);

            if ( dist < (this->simplified_grain_representations[g].getRadius() + userInputs.buffer_between_grains/2.0) ){
                feature_ids[v] = (double)(this->simplified_grain_representations[g].getGrainId());


            }

        }

    }
    
} //endfor ni




scalarvalueType sum_n = constV(0.0);
for (unsigned int i=0; i<userInputs.number_of_variables; i++){
    ni = variable_list.get_scalar_value(i);
    sum_n += ni;
}
for (unsigned int v=0; v<ni.size();v++){
    if (sum_n[v] < 0.01){
        max_op[v] = -1.0;
        feature_ids[v] = -1.0;
    }
}


//n0 = variable_list.get_scalar_value(0);
//double A = this->simplified_grain_representations[0].getGrainId();

//Rho = Rho*n0;
//Rho = std::max(Rho, constV(0.0));

//Rho = A*n0;

// I need to ensure if the value of Rho is negative then the it's set to zero


//std::cout << A << std::endl;

// --- Submitting the terms for the postprocessing expressions ---

pp_variable_list.set_scalar_value_term_RHS(0, feature_ids);
pp_variable_list.set_scalar_value_term_RHS(1, max_op);
//pp_variable_list.set_scalar_value_term_RHS(2, Rho);

}
