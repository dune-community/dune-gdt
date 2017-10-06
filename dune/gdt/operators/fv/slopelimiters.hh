// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_OPERATORS_FV_SLOPELIMITERS_HH
#define DUNE_GDT_OPERATORS_FV_SLOPELIMITERS_HH

#include <dune/xt/common/float_cmp.hh>

namespace Dune {
namespace GDT {


enum class SlopeLimiters
{
  minmod,
  mc,
  superbee,
  no_slope
};


namespace internal {


template <SlopeLimiters slope_limiter>
struct ChooseLimiter
{
  template <class VectorType>
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& slope_center);
};

template <>
struct ChooseLimiter<SlopeLimiters::minmod>
{
  template <class VectorType>
  static VectorType
  limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& /*slope_center*/)
  {
    VectorType ret(0.);
    for (size_t ii = 0; ii < slope_left.size(); ++ii)
      if (slope_left[ii] * slope_right[ii] > 0) // check for equal sign
        ret[ii] = (std::abs(slope_left[ii]) < std::abs(slope_right[ii])) ? slope_left[ii] : slope_right[ii];
    return ret;
  }
};

template <>
struct ChooseLimiter<SlopeLimiters::superbee>
{
  template <class VectorType>
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& slope_center)
  {
    auto slope_left_twice = slope_left;
    slope_left_twice *= 2.0;
    auto slope_right_twice = slope_right;
    slope_right_twice *= 2.0;
    typedef ChooseLimiter<SlopeLimiters::minmod> MinmodType;
    return maxmod(MinmodType::limit(slope_left, slope_right_twice, slope_center),
                  MinmodType::limit(slope_left_twice, slope_right, slope_center));
  }

  template <class VectorType>
  static VectorType maxmod(const VectorType& slope_left, const VectorType& slope_right)
  {
    VectorType ret(0.);
    for (size_t ii = 0; ii < slope_left.size(); ++ii) {
      const auto slope_left_abs = std::abs(slope_left[ii]);
      const auto slope_right_abs = std::abs(slope_right[ii]);
      if (slope_left_abs > slope_right_abs && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_left[ii];
      else if (Dune::XT::Common::FloatCmp::le(slope_left_abs, slope_right_abs) && slope_left[ii] * slope_right[ii] > 0)
        ret[ii] = slope_right[ii];
      else
        ret[ii] = 0.0;
    }
    return ret;
  }
};

template <>
struct ChooseLimiter<SlopeLimiters::mc>
{
  template <class VectorType>
  static VectorType limit(VectorType slope_left, VectorType slope_right, const VectorType& slope_center)
  {
    typedef ChooseLimiter<SlopeLimiters::minmod> MinmodType;
    slope_left *= 2.0;
    slope_right *= 2.0;
    return MinmodType::limit(MinmodType::limit(slope_left, slope_right, slope_center), slope_center, slope_center);
  }
};

template <>
struct ChooseLimiter<SlopeLimiters::no_slope>
{
  template <class VectorType>
  static VectorType
  limit(const VectorType& /*slope_left*/, const VectorType& /*slope_right*/, const VectorType& /*slope_center*/)
  {
    return VectorType(0.);
  }
};


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_SLOPELIMITERS_HH
