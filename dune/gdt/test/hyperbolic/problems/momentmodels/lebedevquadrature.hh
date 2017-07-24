// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/math.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace internal {


struct LebedevHelper
{
  static void get(std::vector<Dune::QuadraturePoint<double, 3>>& quad_vec, const size_t order);
};


} // namespace internal


// The tabulated values are in cartesian coordinates (x, y, z). If cartesian is false, the
// quadrature points are converted to spherical coordinates (\theta, \varphi), with
// 0 \leq \theta \leq \pi and 0 \leq \varphi \leq 2\pi.
template <class FieldType, bool cartesian = true>
class LebedevQuadrature
{
public:
  static Dune::QuadratureRule<FieldType, 2 + cartesian> get(const size_t requested_order)
  {
    size_t index = -1;
    for (size_t ii = 0; ii < allowed_orders_.size(); ++ii) {
      if (allowed_orders_[ii] >= requested_order) {
        index = ii;
        break;
      }
    }
    if (index == size_t(-1))
      std::cerr << "Warning: Requested Lebedev quadrature with order " << requested_order
                << " is not available, using highest available order " << allowed_orders_.back() << "." << std::endl;
    size_t order = (index == size_t(-1)) ? allowed_orders_.back() : allowed_orders_[index];

    std::vector<Dune::QuadraturePoint<double, 3>> quad_vector;
    internal::LebedevHelper::get(quad_vector, order);

    Dune::QuadratureRule<FieldType, 3> quad_rule;
    for (const auto& quad_point : quad_vector)
      quad_rule.emplace_back(quad_point);

    return helper<>::create(quad_rule);
  }

private:
  template <bool cart = cartesian, class anything = void>
  struct helper
  {
    static Dune::QuadratureRule<FieldType, 3> create(const Dune::QuadratureRule<FieldType, 3> quad_in)
    {
      return quad_in;
    }
  };

  template <class anything>
  struct helper<false, anything>
  {
    static Dune::QuadratureRule<FieldType, 2> create(const Dune::QuadratureRule<FieldType, 3> quad_in)
    {
      Dune::QuadratureRule<FieldType, 3> ret;
      for (const auto& point : quad_in) {
        ret.emplace_back(Dune::QuadraturePoint<FieldType, 3>(
            XT::Common::CoordinateConverter<FieldType>::to_spherical(point.position()), point.weight()));
      }
      return ret;
    }
  };

  static const std::vector<size_t> allowed_orders_;
}; // class LebedevQuadrature

template <class FieldType, bool cartesian>
const std::vector<size_t> LebedevQuadrature<FieldType, cartesian>::allowed_orders_ = {
    3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
    41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH
