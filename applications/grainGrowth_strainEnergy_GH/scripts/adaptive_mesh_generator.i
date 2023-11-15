/*

adaptive_list_generator.i

//Works for periodic boundary conditions ONLY
Generates a list of subsections for the adaptivity subsection for
"n"
*/

ngrains=38;

fl=create("adaptive_mesh_text.dat");
for (j=0;j<ngrains;j++){
  write, fl, format="%s", "subsection Refinement criterion: n";
  write, fl, format="%i\n", j;
  write, fl, format="%s\n", "# Select whether the mesh is refined based on the variable value (VALUE),";
  write, fl, format="%s\n", "# its gradient (GRADIENT), or both (VALUE_AND_GRADIENT)";
  write, fl, format="%s\n", "set Criterion type = VALUE";
  write, fl, format="%s\n", "# Set the lower and upper bounds for the value-based refinement window";
  write, fl, format="%s\n", "set Value lower bound = 0.001";
  write, fl, format="%s\n", "set Value upper bound = 0.999";
  write, fl, format="%s\n", "end";
  write, fl, format="%s\n", "";
}
close, fl;
