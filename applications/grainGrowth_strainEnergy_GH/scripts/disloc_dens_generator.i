/*

ic_generator.i

Calculates a vector of normally distributed dislocation density values

*/

//Calculating rescaled dislocation density values from initial conditions for PRL paper

N=38;

rho_ind=array(0,N);
rho_units=array(0.0,N);
f1=open("DislocationDensity.txt");
for(j=1;j<=N;j++)
  read, f1, rho_ind(j), rho_units(j);
close, f1;

rho_mean=avg(rho_units);

//Rescaled dislocation density
rho_i=rho_units/rho_mean;

//Writing the vector for densities

fl=create("dislocation_densities.txt");

write, fl, format="%s", "{";
for(j=1;j<=N;j++){
    if (j<N){
        if (j%4) {
            write, fl, format="%.5f, ", rho_i(j);
        } else {
            write, fl, format="%.5f, \n", rho_i(j);
        }
    }
    if (j==N) write, fl, format="%.5f", rho_i(j);
}
write, fl, format="%s", "};";

close, fl;
