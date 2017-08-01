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


template<size_t order>
struct LebedevHelper
{
  static void get(std::vector<Dune::QuadraturePoint<double, 3>>& quad_vec);
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
    switch (order) {
      case 3: internal::LebedevHelper<3>::get(quad_vector); break;
      case 5: internal::LebedevHelper<5>::get(quad_vector); break;
      case 7: internal::LebedevHelper<7>::get(quad_vector); break;
      case 9: internal::LebedevHelper<9>::get(quad_vector); break;
      case 11: internal::LebedevHelper<11>::get(quad_vector); break;
      case 13: internal::LebedevHelper<13>::get(quad_vector); break;
      case 15: internal::LebedevHelper<15>::get(quad_vector); break;
      case 17: internal::LebedevHelper<17>::get(quad_vector); break;
      case 19: internal::LebedevHelper<19>::get(quad_vector); break;
      case 21: internal::LebedevHelper<21>::get(quad_vector); break;
      case 23: internal::LebedevHelper<23>::get(quad_vector); break;
      case 25: internal::LebedevHelper<25>::get(quad_vector); break;
      case 27: internal::LebedevHelper<27>::get(quad_vector); break;
      case 29: internal::LebedevHelper<29>::get(quad_vector); break;
      case 31: internal::LebedevHelper<31>::get(quad_vector); break;
      case 35: internal::LebedevHelper<35>::get(quad_vector); break;
      case 41: internal::LebedevHelper<41>::get(quad_vector); break;
      case 47: internal::LebedevHelper<47>::get(quad_vector); break;
      case 53: internal::LebedevHelper<53>::get(quad_vector); break;
      case 59: internal::LebedevHelper<59>::get(quad_vector); break;
      case 65: internal::LebedevHelper<65>::get(quad_vector); break;
      case 71: internal::LebedevHelper<71>::get(quad_vector); break;
      case 77: internal::LebedevHelper<77>::get(quad_vector); break;
      case 83: internal::LebedevHelper<83>::get(quad_vector); break;
      case 89: internal::LebedevHelper<89>::get(quad_vector); break;
      case 95: internal::LebedevHelper<95>::get(quad_vector); break;
      case 101: internal::LebedevHelper<101>::get(quad_vector); break;
      case 107: internal::LebedevHelper<107>::get(quad_vector); break;
      case 113: internal::LebedevHelper<113>::get(quad_vector); break;
      case 119: internal::LebedevHelper<119>::get(quad_vector); break;
      case 125: internal::LebedevHelper<125>::get(quad_vector); break;
      case 131: internal::LebedevHelper<131>::get(quad_vector); break;
    default: DUNE_THROW(NotImplemented, "Requested order is not available!");
    }

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
