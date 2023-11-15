/*

dependency_list_generator.i

//Works for periodic boundary conditions ONLY
Generates a list of order parameters separated by commas and quotations signs
"n"
*/


//Number of order parameters
nops = 38;

//Writing out positions readable in icc
fl=create("list_of_ops.txt");
write, fl, "List of order parameters";
write, fl, format="%s", "\""
for(j=0;j<nops;j++){
  write, fl, format="%s", "n";
  write, fl, format="%d", j;
  //if ((j+1)%4 == 0){
    if (j==nops-1){
      write, fl, format="%s", "\"";
    } else {
      write, fl, format="%s", ", ";
    }
  //} else {
  //    write, fl, format="%s", ", ";
  //}
}
close, fl;
