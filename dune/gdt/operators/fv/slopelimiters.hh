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

#include <dune/xt/la/container/eigen.hh>

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


template <SlopeLimiters slope_limiter, class VectorType>
struct ChooseLimiter
{
  static VectorType limit(const VectorType& u_entity, const VectorType& u_left, const VectorType& u_right);
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::minmod, VectorType>
{
  typedef typename XT::Common::VectorAbstraction<VectorType> Abstr;
  static VectorType limit(const VectorType& u_entity, const VectorType& u_left, const VectorType& u_right)
  {
    VectorType ret(u_entity.size(), 0., 0);
    for (size_t ii = 0; ii < u_entity.size(); ++ii) {
      const auto slope_left = Abstr::get_entry(u_entity, ii) - Abstr::get_entry(u_left, ii);
      const auto slope_right = Abstr::get_entry(u_right, ii) - Abstr::get_entry(u_entity, ii);
      const auto slope_center = 0.5 * (Abstr::get_entry(u_right, ii) - Abstr::get_entry(u_left, ii));
      // check for equal sign
      if (slope_left * slope_right > 0 && slope_center * slope_right > 0) {
        const auto slope_left_abs = std::abs(slope_left);
        const auto slope_right_abs = std::abs(slope_right);
        const auto slope_center_abs = std::abs(slope_center);
        if (XT::Common::FloatCmp::lt(slope_left_abs, slope_right_abs)) {
          if (XT::Common::FloatCmp::lt(slope_left_abs, slope_center_abs))
            Abstr::set_entry(ret, ii, slope_left);
        } else if (XT::Common::FloatCmp::lt(slope_right_abs, slope_center_abs))
          Abstr::set_entry(ret, ii, slope_right);
        else
          Abstr::set_entry(ret, ii, slope_center);
      }
    }
    return ret;
  }
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::superbee, VectorType>
{
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope)
  {
    typedef ChooseLimiter<SlopeLimiters::minmod, VectorType> MinmodType;
    return maxmod(MinmodType::limit(slope_left, slope_right * 2.0, centered_slope),
                  MinmodType::limit(slope_left * 2.0, slope_right, centered_slope));
  }

  static VectorType maxmod(const VectorType& slope_left, const VectorType& slope_right)
  {
    VectorType ret;
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

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::mc, VectorType>
{
  static VectorType limit(const VectorType& slope_left, const VectorType& slope_right, const VectorType& centered_slope)
  {
    typedef ChooseLimiter<SlopeLimiters::minmod, VectorType> MinmodType;
    return MinmodType::limit(
        MinmodType::limit(slope_left * 2.0, slope_right * 2.0, centered_slope), centered_slope, centered_slope);
  }
};

template <class VectorType>
struct ChooseLimiter<SlopeLimiters::no_slope, VectorType>
{
  typedef typename XT::Common::VectorAbstraction<VectorType> Abstr;
  static VectorType limit(const VectorType& u_entity, const VectorType& /*u_left*/, const VectorType& /*u_right*/)
  {
    return VectorType(u_entity.size(), 0., 0);
  }
};


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_SLOPELIMITERS_HH
