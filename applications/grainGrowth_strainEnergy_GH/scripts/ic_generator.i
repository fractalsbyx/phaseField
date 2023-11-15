/*

ic_generator.i

//Works for periodic boundary conditions ONLY
generates an array of random positions and radii of nuclei in 2D or 3D
in the following format
Positions:
 {{0.1,0.3,0},{0.8,0.7,0},{0.5,0.2,0},{0.4,0.4,0},{0.3,0.9,0},...,{1,1,0},{0.7,0.95,0}};
 Radii:
 {12, 14, 19, 16, 11, 12, 17, ... , 11, 14};

*/

//Calculating random seed positions from initial conditions for PRL paper

//MAIN
//Defining constants

//dimensionality
dim=2;

//System size
lx=6.305;
ly=6.305;

lvec=[lx,ly];

//Number of seeds
nseeds = 10;
//delta x
delx=0.025;
//Average radius of seeds
nrad = 3.5*delx;

//Initializing position and dimension vectors;
nucposvec_inds=array(0,[2,dim,nseeds]);
nucposvec=array(0.0,[2,dim,nseeds]);

f1=open("rand_locs_inds.txt");
read, f1, nucposvec_inds;
close, f1;

//nucposvec(1,)=1.0*nucposvec_inds(2,)*delx;
//nucposvec(2,)=ly-1.0*nucposvec_inds(1,)*delx;

nucposvec(1,)=1.0*nucposvec_inds(1,)*delx;
nucposvec(2,)=1.0*nucposvec_inds(2,)*delx;


if (dim==2){
    winkill, 0;
    window, 0, dpi=100;
    limits, 0, lx, 0, ly;
    fma;
    for(j=1;j<=nseeds;j++){
        plmk, nucposvec(2,j),nucposvec(1,j), marker=4, msize=0.5, width=10;
        //rdline;
    }
}

//Writing out positions readable in icc
fl=create("nuclei_data.txt");
write, fl, "Nuclei coordinates";
write, fl, format="%s", "{";
for(j=1;j<=nseeds;j++){
    if (j<nseeds){
        if (j%4) {
            write, fl, format="%s", "{";
            for (k=1;k<=dim;k++){
                if (k<dim) write, fl, format="%.3e, ", nucposvec(k,j);
                if (k==dim) write, fl, format="%.3e", nucposvec(k,j);
            }
            write, fl, format="%s", "},";
        } else {
            write, fl, format="%s", "{";
            for (k=1;k<=dim;k++){
                if (k<dim) write, fl, format="%.3e, ", nucposvec(k,j);
                if (k==dim) write, fl, format="%.3e", nucposvec(k,j);
            }
            write, fl, format="%s\n", "},";
        }
    }
    if (j==nseeds){
        write, fl, format="%s", "{";
        for (k=1;k<=dim;k++){
            if (k<dim) write, fl, format="%.3e, ", nucposvec(k,j);
            if (k==dim) write, fl, format="%.3e", nucposvec(k,j);
        }
        write, fl, format="%s", "}";
    }
}
write, fl, format="%s \n\n", "};";
