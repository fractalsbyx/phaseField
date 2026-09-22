// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <cmath>

/**
 * This file provides some operations on VectorizedArray that are not provided by
 * dealii. Use the std namespace so we can call them with std::function
 */
namespace std
{
  // NOLINTBEGIN(cert-dcl58-cpp, readability-identifier-naming,
  // readability-identifier-length)

  template <typename Number, std::size_t width>
  inline ::VectorizedArray<Number, width>
  erf(const ::VectorizedArray<Number, width> &x)
  {
    ::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::erf(x[i]);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::VectorizedArray<Number, width>
  erfc(const ::VectorizedArray<Number, width> &x)
  {
    ::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::erfc(x[i]);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::VectorizedArray<Number, width>
  atan2(const ::VectorizedArray<Number, width> &y,
        const ::VectorizedArray<Number, width> &x)
  {
    ::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::atan2(y[i], x[i]);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::VectorizedArray<Number, width>
  fmod(const ::VectorizedArray<Number, width> &numer, const Number denom)
  {
    ::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::fmod(numer[i], denom);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::VectorizedArray<Number, width>
  fmod(const ::VectorizedArray<Number, width> &numer,
       const ::VectorizedArray<Number, width> &denom)
  {
    ::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::fmod(numer[i], denom[i]);
    return out;
  }

  // NOLINTEND(cert-dcl58-cpp, readability-identifier-naming,
  // readability-identifier-length)

} // namespace std
