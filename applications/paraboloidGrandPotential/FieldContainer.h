#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

template <unsigned int dim>
struct FieldContainer
{
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a);

  scalarValue val = constV(0.);
  scalarGrad  grad;

  FieldContainer<dim>
  operator+(const FieldContainer<dim> &other) const
  {
    return {val + other.val, grad + other.grad};
  }

  FieldContainer<dim>
  operator-(const FieldContainer<dim> &other) const
  {
    return {val - other.val, grad - other.grad};
  }

  FieldContainer<dim> &
  operator+=(const FieldContainer<dim> &other)
  {
    val += other.val;
    grad += other.grad;
    return *this;
  }

  FieldContainer<dim> &
  operator-=(const FieldContainer<dim> &other)
  {
    val -= other.val;
    grad -= other.grad;
    return *this;
  }

  FieldContainer<dim>
  operator+() const
  {
    return *this;
  }

  FieldContainer<dim>
  operator-() const
  {
    return {-val, -grad};
  }

  FieldContainer<dim>
  operator*(const scalarValue &scalar) const
  {
    return {val * scalar, grad * scalar};
  }

  FieldContainer<dim>
  operator*(const double &scalar) const
  {
    return {val * scalar, grad * scalar};
  }

  FieldContainer<dim>
  operator*(const FieldContainer<dim> &other) const
  {
    return {val * other.val, grad * other.val + val * other.grad};
  }

  FieldContainer<dim> &
  operator*=(const scalarValue &scalar)
  {
    val *= scalar;
    grad *= scalar;
    return *this;
  }

  FieldContainer<dim> &
  operator*=(const double &scalar)
  {
    val *= scalar;
    grad *= scalar;
    return *this;
  }

  FieldContainer<dim> &
  operator*=(const FieldContainer<dim> &other)
  {
    grad = grad * other.val + val * other.grad;
    val *= other.val;
    return *this;
  }

  FieldContainer<dim>
  operator/(const scalarValue &scalar) const
  {
    return {val / scalar, grad / scalar};
  }

  FieldContainer<dim>
  operator/(const double &scalar) const
  {
    return {val / scalar, grad / scalar};
  }

  FieldContainer<dim>
  operator/(const FieldContainer<dim> &other) const
  {
    return {val / other.val,
            (grad * other.val - val * other.grad) / (other.val * other.val)};
  }

  FieldContainer<dim> &
  operator/=(const scalarValue &scalar)
  {
    val /= scalar;
    grad /= scalar;
    return *this;
  }

  FieldContainer<dim> &
  operator/=(const double &scalar)
  {
    val /= scalar;
    grad /= scalar;
    return *this;
  }

  FieldContainer<dim> &
  operator/=(const FieldContainer<dim> &other)
  {
    val /= other.val;
    grad = (grad * other.val - val * other.grad) / (other.val * other.val);
    return *this;
  }

  static FieldContainer<dim>
  field_x_variation(const FieldContainer<dim> &field,
                    const FieldContainer<dim> &variation)
  {
    return {field.val * variation.val - field.grad * variation.grad,
            field.val * variation.grad};
  }
};