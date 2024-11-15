#ifndef FIELDCONTAINER_H
#define FIELDCONTAINER_H

#include <deal.II/base/vectorization.h>

#include "../../include/variableContainer.h"

template <int dim>
struct FieldContainer
{
  dealii::VectorizedArray<double> val = dealii::make_vectorized_array(0.);
  dealii::Tensor<1, dim, dealii::VectorizedArray<double>> grad;

  inline FieldContainer<dim>
  operator*(const double &other) const
  {
    return FieldContainer<dim> {other * val, other * grad};
  }

  // Use with caution
  inline FieldContainer<dim>
  operator+(const FieldContainer<dim> &other) const
  {
    return FieldContainer<dim> {val + other.val, grad + other.grad};
  }

  // Use with caution
  inline FieldContainer<dim>
  operator-(const FieldContainer<dim> &other) const
  {
    return FieldContainer<dim> {val - other.val, grad - other.grad};
  }

  inline FieldContainer<dim>
  operator+() const
  {
    return FieldContainer<dim> {val, grad};
  }

  inline FieldContainer<dim>
  operator-() const
  {
    return FieldContainer<dim> {-val, -grad};
  }

  inline void
  operator+=(const FieldContainer<dim> &other)
  {
    val += other.val;
    grad += other.grad;
  }

  inline void
  operator-=(const FieldContainer<dim> &other)
  {
    val -= other.val;
    grad -= other.grad;
  }

  inline void
  operator*=(const double &other)
  {
    val *= other;
    grad *= other;
  }

  inline FieldContainer<dim>
  inverse() const
  {
    return FieldContainer<dim> {dealii::make_vectorized_array(1.) / val,
                                -grad / (val * val)};
  }
};

template <int dim>
inline FieldContainer<dim>
field_x_field(const FieldContainer<dim> &field1, const FieldContainer<dim> &field2)
{
  return FieldContainer<dim> {field1.val * field2.val,
                              field1.grad * field2.val + field1.val * field2.grad};
}

template <int dim>
inline FieldContainer<dim>
field_x_variation(const FieldContainer<dim> &field, const FieldContainer<dim> &variation)
{
  return FieldContainer<dim> {field.val * variation.val - field.grad * variation.grad,
                              field.val * variation.grad};
}

template <int dim>
inline FieldContainer<dim>
operator*(const double &val, const FieldContainer<dim> &field)
{
  return FieldContainer<dim> {val * field.val, val * field.grad};
}

#endif