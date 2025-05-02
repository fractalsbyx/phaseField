#ifndef LINALG_H
#define LINALG_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "FieldContainer.h"

template <typename number>
using boost_vector = boost::numeric::ublas::vector<number>;
template <typename number>
using boost_symmet = boost::numeric::ublas::symmetric_matrix<number>;
template <typename number>
using boost_matrix = boost::numeric::ublas::matrix<number>;

// Returns true on success, false if matrix is not SPD
// Step 1: Cholesky decomposition
template <typename number>
bool
cholesky_decompose(const boost_symmet<number> &A, boost_matrix<number> &L)
{
  using std::sqrt;
  const std::size_t n = A.size1();

  if (A.size1() != A.size2())
    {
      throw std::invalid_argument("Matrix must be square");
    }

  L = boost_matrix<number>(n, n);

  for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j <= i; ++j)
        {
          number sum = A(i, j);
          for (std::size_t k = 0; k < j; ++k)
            {
              sum -= L(i, k) * L(j, k);
            }

          if (i == j)
            {
              // if (sum <= 0.0)
              //   return false; // Not positive definite
              L(i, j) = sqrt(sum);
            }
          else
            {
              L(i, j) = sum / L(j, j);
            }
        }
    }
  return true;
}

// Step 2: Forward substitution
template <typename number>
boost_vector<number>
forward_substitution(const boost_matrix<number> &L, const boost_vector<number> &b)
{
  std::size_t          n = L.size1();
  boost_vector<number> y(n);
  for (std::size_t i = 0; i < n; ++i)
    {
      number sum = b(i);
      for (std::size_t j = 0; j < i; ++j)
        {
          sum -= L(i, j) * y(j);
        }
      y(i) = sum / L(i, i);
    }
  return y;
}

// Step 3: Backward substitution (for L^T)
template <typename number>
boost_vector<number>
backward_substitution(const boost_matrix<number> &L, const boost_vector<number> &y)
{
  int                  n = L.size1();
  boost_vector<number> x(n);
  for (int i = n - 1; i >= 0; --i)
    {
      number sum = y(i);
      for (int j = i + 1; j < n; ++j)
        {
          sum -= L(j, i) * x(j); // Note: L^T used here
        }
      x(i) = sum / L(i, i);
    }
  return x;
}

// Step 4: Compute A inverse
template <typename number>
boost_matrix<number>
inverse_from_cholesky(const boost_symmet<number> &A)
{
  const std::size_t    n = A.size1();
  boost_matrix<number> L;
  if (!cholesky_decompose(A, L))
    {
      throw std::runtime_error("Matrix is not positive definite");
    }

  boost_matrix<number> Ainv(n, n);
  for (std::size_t i = 0; i < n; ++i)
    {
      // i-th column of identity
      boost_vector<number> e(n);
      e(i) = number(1.0);

      // Solve L y = e
      boost_vector<number> y = forward_substitution(L, e);
      // Solve L^T x = y
      boost_vector<number> x = backward_substitution(L, y);

      // Set i-th column of Ainv
      for (std::size_t j = 0; j < n; ++j)
        {
          Ainv(j, i) = x(j);
        }
    }

  return Ainv;
}

#endif // LINALG_H