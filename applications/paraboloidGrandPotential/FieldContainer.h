#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

/**
 * @brief Class for holding fields and their spatial derivatives OR variations in their
 * strong form
 * @details Variations are in strong form: [val + div(grad)]
 * @tparam dim The dimension of the problem
 */
template <unsigned int dim>
struct FieldContainer
{
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)

  scalarValue val  = constV(0.);
  scalarGrad  grad = {};

  FieldContainer<dim>
  operator+(const FieldContainer<dim> &other) const
  {
    return {val + other.val, grad + other.grad};
  }

  FieldContainer<dim>
  operator+(const scalarValue &constant) const
  {
    return {val + constant, grad};
  }

  FieldContainer<dim>
  operator+(const double &constant) const
  {
    return {val + constant, grad};
  }

  FieldContainer<dim>
  operator+() const
  {
    return *this;
  }

  template <typename number>
  FieldContainer<dim> &
  operator+=(const number &other)
  {
    *this = *this + other;
    return *this;
  }

  FieldContainer<dim>
  operator-(const FieldContainer<dim> &other) const
  {
    return {val - other.val, grad - other.grad};
  }

  FieldContainer<dim>
  operator-() const
  {
    return {-val, -grad};
  }

  FieldContainer<dim>
  operator-(const scalarValue &constant) const
  {
    return {val - constant, grad};
  }

  FieldContainer<dim>
  operator-(const double &constant) const
  {
    return {val - constant, grad};
  }

  template <typename number>
  FieldContainer<dim> &
  operator-=(const number &other)
  {
    *this = *this - other;
    return *this;
  }

  FieldContainer<dim>
  operator*(const scalarValue &constant) const
  {
    return {val * constant, grad * constant};
  }

  FieldContainer<dim>
  operator*(const double &constant) const
  {
    return {val * constant, grad * constant};
  }

  FieldContainer<dim>
  operator*(const FieldContainer<dim> &other) const
  {
    return {val * other.val, val * other.grad + grad * other.val};
  }

  template <typename number>
  FieldContainer<dim> &
  operator*=(const number &other)
  {
    *this = *this * other;
    return *this;
  }

  FieldContainer<dim>
  operator/(const scalarValue &constant) const
  {
    return {val / constant, grad / constant};
  }

  FieldContainer<dim>
  operator/(const double &constant) const
  {
    return {val / constant, grad / constant};
  }

  FieldContainer<dim>
  operator/(const FieldContainer<dim> &other) const
  {
    return {val / other.val,
            (grad * other.val - val * other.grad) / (other.val * other.val)};
  }

  template <typename number>
  FieldContainer<dim> &
  operator/=(const number &other)
  {
    *this = *this / other;
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

template <typename number, unsigned int dim>
FieldContainer<dim>
operator+(const number &other, const FieldContainer<dim> &field)
{
  return field + other;
}

template <typename number, unsigned int dim>
FieldContainer<dim>
operator-(const number &other, const FieldContainer<dim> &field)
{
  return -field + other;
}

template <typename number, unsigned int dim>
FieldContainer<dim>
operator*(const number &other, const FieldContainer<dim> &field)
{
  return field * other;
}

template <typename number, unsigned int dim>
FieldContainer<dim>
operator/(const number &other, const FieldContainer<dim> &field)
{
  FieldContainer<dim> one = {dealii::make_vectorized_array(1.0), {}};
  return (one / field) * other;
}

template <unsigned int dim>
struct Variation
{
  using scalarValue = dealii::VectorizedArray<double>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<double>>;
#define constV(a) dealii::make_vectorized_array(a)

  scalarValue val = constV(0.);
  scalarGrad  vec = {};

  Variation<dim>
  operator+(const Variation<dim> &other) const
  {
    return {val + other.val, vec + other.vec};
  }

  Variation<dim>
  operator+(const scalarValue &scalar) const
  {
    return {val + scalar, vec};
  }

  Variation<dim>
  operator+(const double &scalar) const
  {
    return {val + scalar, vec};
  }

  Variation<dim>
  operator+() const
  {
    return *this;
  }

  template <typename number>
  Variation<dim> &
  operator+=(const number &other)
  {
    *this = *this + other;
    return *this;
  }

  Variation<dim>
  operator-(const Variation<dim> &other) const
  {
    return {val - other.val, vec - other.vec};
  }

  Variation<dim>
  operator-() const
  {
    return {-val, -vec};
  }

  Variation<dim>
  operator-(const scalarValue &scalar) const
  {
    return {val - scalar, vec};
  }

  Variation<dim>
  operator-(const double &scalar) const
  {
    return {val - scalar, vec};
  }

  template <typename number>
  Variation<dim> &
  operator-=(const number &other)
  {
    *this = *this - other;
    return *this;
  }

  Variation<dim>
  operator*(const scalarValue &constant_scalar) const
  {
    return {val * constant_scalar, vec * constant_scalar};
  }

  Variation<dim>
  operator*(const double &constant_scalar) const
  {
    return {val * constant_scalar, vec * constant_scalar};
  }

  template <typename number>
  Variation<dim> &
  operator*=(const number &other)
  {
    *this = *this * other;
    return *this;
  }

  Variation<dim>
  operator/(const scalarValue &constant_scalar) const
  {
    return {val / constant_scalar, vec / constant_scalar};
  }

  Variation<dim>
  operator/(const double &constant_scalar) const
  {
    return {val / constant_scalar, vec / constant_scalar};
  }

  template <typename number>
  Variation<dim> &
  operator/=(const number &other)
  {
    *this = *this / other;
    return *this;
  }
};

template <typename number, unsigned int dim>
Variation<dim>
operator+(const number &other, const Variation<dim> &variation)
{
  return variation + other;
}

template <typename number, unsigned int dim>
Variation<dim>
operator-(const number &other, const Variation<dim> &variation)
{
  return -variation + other;
}

template <typename number, unsigned int dim>
Variation<dim>
operator*(const number &other, const Variation<dim> &variation)
{
  return variation * other;
}

/**
 * @brief Multiply a field and a variation
 * @param field The field to multiply
 * @param variation The variation to multiply
 */
template <unsigned int dim>
Variation<dim>
operator*(const FieldContainer<dim> &field, const Variation<dim> &variation)
{
  return {field.val * variation.val - field.grad * variation.vec,
          field.val * variation.vec};
}

/**
 * @brief Multiply a field and a variation
 * @param field The field to multiply
 * @param variation The variation to multiply
 */
template <unsigned int dim>
Variation<dim>
operator*(const Variation<dim> &variation, const FieldContainer<dim> &field)
{
  return field * variation;
}