// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/config.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <prismspf/utilities/vectorized_operations.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

namespace Utilities::SignedDistance
{

  /**
   * Compute the signed distance from a point to a sphere.
   *
   * @param p The point to compute the signed distance from.
   * @param center The center of the sphere.
   * @param radius The radius of the sphere.
   * @return The signed distance from the point to the sphere.
   */
  template <int dim>
  inline double
  sdf_sphere(const dealii::Point<dim> &p, const dealii::Point<dim> &center, double radius)
  {
    return (p - center).norm() - radius;
  }

  /**
   */
  template <int dim>
  inline double
  qsdf_sphere(const dealii::Point<dim> &p,
              const dealii::Point<dim> &center,
              double                    radius)
  {
    return 0.5 * radius * radius - (p - center).norm_square() / radius;
  }

  /**
   */
  template <int dim>
  inline double
  sdf_plane(const dealii::Point<dim>     &p,
            const dealii::Point<dim>     &plane_point,
            const dealii::Tensor<1, dim> &normal)
  {
    return (p - plane_point) * normal / normal.norm();
  }

  /**
   */
  inline double
  sdf_cylinder(const dealii::Point<3>     &p,
               const dealii::Point<3>     &axis_point,
               const dealii::Tensor<1, 3> &axis_direction,
               double                      radius)
  {
    dealii::Tensor<1, 3> r1 = p - axis_point;
    dealii::Tensor<1, 3> parallel =
      (r1 * axis_direction / axis_direction.norm()) * axis_direction;
    return (r1 - parallel).norm() - radius;
  }
} // namespace Utilities::SignedDistance

PRISMS_PF_END_NAMESPACE
