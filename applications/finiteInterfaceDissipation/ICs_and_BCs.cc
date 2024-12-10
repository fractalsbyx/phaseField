// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  // Precalculating everything makes writing initial conditions easier. May take slightly
  // more runtime.
  // Custom coordinate system
  double center[3] = {0.5 * userInputs.domain_size[0],
                      0.5 * userInputs.domain_size[1],
                      (dim > 2) * userInputs.domain_size[2]};
  double x         = p[0] - center[0];
  double y         = p[1] - center[1];
  double z         = (dim < 3) ? 0.0 : p[2] - center[2];
  double r2        = x * x + y * y + z * z;

  // TODO: make order parameters
  double p1 = interface(0.75 * (0.5 * (r0 * r0 - x * x) / r0) / Sys.eta);
  double p2 = 1.0 - p1;
  scalar_IC = 0.5;
  if (index == 0)
    {
      scalar_IC = p1;
    }
  if (index == 1)
    {
      scalar_IC = 0.10;
    }
  if (index == 2)
    {
      scalar_IC = 0.30;
    }
  if (index == 3)
    {
      scalar_IC = p2;
    }
  if (index == 4)
    {
      scalar_IC = 0.90;
    }
  if (index == 5)
    {
      scalar_IC = 0.05;
    }

  // ===========================================================================
  // Submit fields
  // ===========================================================================
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
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
