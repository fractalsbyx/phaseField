// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

#include <array>
#include <numbers>

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // ---------------------------------------------------------------------
  std::map<uint, uint>        op_index;
  std::map<std::string, uint> comp_index;
  uint                        var_index = 0;
  for (const auto &comp_name : isoSys.comp_names)
    {
      comp_index[comp_name] = var_index++;
    }
  uint eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : isoSys.order_params)
    {
      op_index[eta_index++] = var_index++;
    }

  // Custom coordinate system
  double center[3] = {0.5 * userInputs.domain_size[0],
                      0.5 * userInputs.domain_size[1],
                      (dim > 2) * userInputs.domain_size[2]};
  double x         = p[0] - center[0];
  double y         = p[1] - center[1];
  double z         = (dim < 3) ? 0.0 : p[2] - center[2];
  double r2        = x * x + y * y + z * z;
  (void) r2;

  // ---------------------------------------------------------------------
  // TODO: ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------

  // Make relevant geometries
  [[maybe_unused]] double circular     = interface(0.5 * (r0 * r0 - r2) / r0);
  [[maybe_unused]] double flat         = interface(0.5 * (r0 * r0 - y * y) / r0);
  [[maybe_unused]] double bottom_strip = interface(r0 - p[2]);

  constexpr std::array<double, 3> a = {
    {0.0, 1.0, 2.0}
  };
  constexpr double                     tau    = std::acos(-1.0) * 2.0;
  double                               s      = 2.0 / (3.0 * std::sqrt(3.0));
  double                               offset = 1.0;
  double                               m      = 0.0;
  std::array<std::array<double, 2>, 3> b      = {
    {{{std::cos(tau * a[0] / 3.0), std::sin(tau * a[0] / 3.0)}},
     {{std::cos(tau * a[1] / 3.0), std::sin(tau * a[1] / 3.0)}},
     {{std::cos(tau * a[2] / 3.0), std::sin(tau * a[2] / 3.0)}}}
  };
  std::array<std::array<double, 2>, 3> d = {
    {
     {{std::cos(tau * a[0] / 3.0 + tau / 12.0),
        std::sin(tau * a[0] / 3.0 + tau / 12.0)}},
     {{std::cos(tau * a[1] / 3.0 + tau / 12.0),
        std::sin(tau * a[1] / 3.0 + tau / 12.0)}},
     {{std::cos(tau * a[2] / 3.0 + tau / 12.0),
        std::sin(tau * a[2] / 3.0 + tau / 12.0)}} //
    }
  };

  std::array<double, 3> w = {
    {0.0, 0.0, 0.0}
  };
  for (uint i = 0; i < 3; i++)
    {
      for (uint j = 0; j < 3; j++)
        {
          w[i] += std::cos(tau * ((x / s0 + offset * s * d[i][0]) * b[j][0] +
                                  (y / s0 + offset * s * d[i][1]) * b[j][1] + m));
        }
      w[i] *= s0 / tau;
    }

  // TODO: Populate eta0 with the initial condition for each order parameters
  std::vector<double> eta0(isoSys.order_params.size(), 0.0);
  eta0[0] = 1.0 - bottom_strip;
  // eta0[1] = bottom_strip * interface(eutectic_contour_2D(x, 3, 0, s0));
  // eta0[2] = bottom_strip * interface(eutectic_contour_2D(x, 3, 1, s0));
  // eta0[3] = bottom_strip * interface(eutectic_contour_2D(x, 3, 2, s0));

  // eta0[1] = bottom_strip * ((1.0 / 3.0) + 0.1 * dist(rng));
  // eta0[2] = bottom_strip * ((1.0 / 3.0) + 0.1 * dist(rng));
  // eta0[3] = bottom_strip * ((1.0 / 3.0) + 0.1 * dist(rng));

  eta0[1] = bottom_strip * interface(w[0]);
  eta0[2] = bottom_strip * interface(w[1]);
  eta0[3] = bottom_strip * interface(w[2]);

  // ---------------------------------------------------------------------
  //
  // ---------------------------------------------------------------------

  // Submit the fields
  var_index = 0;
  for (uint comp_index = 0; comp_index < isoSys.comp_names.size(); comp_index++)
    {
      if (index == var_index)
        {
          double mu0 = 0.;
          eta_index  = 0;
          for (const auto &phase_index : isoSys.order_params)
            {
              auto &phase_comp_info = isoSys.phases.at(phase_index).comps.at(comp_index);
              mu0 += eta0[eta_index] * isoSys.Vm * phase_comp_info.k_well *
                     (phase_comp_info.x0 - phase_comp_info.c_min);
              eta_index++;
            }
          scalar_IC = mu0;
        }
      var_index++;
    }
  eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : isoSys.order_params)
    {
      if (index == op_index[eta_index])
        {
          scalar_IC = eta0[eta_index];
        }
      eta_index++;
      var_index++;
    }

  // ---------------------------------------------------------------------
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

  // --------------------------------------------------------------------------
}
