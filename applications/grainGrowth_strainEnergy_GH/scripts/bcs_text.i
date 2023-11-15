/*

BCs_text.i

*/

ngrains=38;

fl=create("bcs_text.dat");
for (j=0;j<ngrains;j++){
  write, fl, format="%s", "set Boundary condition for variable n";
  write, fl, format="%i", j
  write, fl, format="%s\n",  " = NATURAL"
}
close, fl;
