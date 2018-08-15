// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/coordinates.hh>


#include <dune/xt/data/lebedev/lebedev_data.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

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

    std::vector<std::pair<std::array<double, 3>, double>> quad_vector;
    switch (order) {
      case 3:
        internal::LebedevData<3>::get(quad_vector);
        break;
      case 5:
        internal::LebedevData<5>::get(quad_vector);
        break;
      case 7:
        internal::LebedevData<7>::get(quad_vector);
        break;
      case 9:
        internal::LebedevData<9>::get(quad_vector);
        break;
      case 11:
        internal::LebedevData<11>::get(quad_vector);
        break;
      case 13:
        internal::LebedevData<13>::get(quad_vector);
        break;
      case 15:
        internal::LebedevData<15>::get(quad_vector);
        break;
      case 17:
        internal::LebedevData<17>::get(quad_vector);
        break;
      case 19:
        internal::LebedevData<19>::get(quad_vector);
        break;
      case 21:
        internal::LebedevData<21>::get(quad_vector);
        break;
      case 23:
        internal::LebedevData<23>::get(quad_vector);
        break;
      case 25:
        internal::LebedevData<25>::get(quad_vector);
        break;
      case 27:
        internal::LebedevData<27>::get(quad_vector);
        break;
      case 29:
        internal::LebedevData<29>::get(quad_vector);
        break;
      case 31:
        internal::LebedevData<31>::get(quad_vector);
        break;
      case 35:
        internal::LebedevData<35>::get(quad_vector);
        break;
      case 41:
        internal::LebedevData<41>::get(quad_vector);
        break;
      case 47:
        internal::LebedevData<47>::get(quad_vector);
        break;
      case 53:
        internal::LebedevData<53>::get(quad_vector);
        break;
      case 59:
        internal::LebedevData<59>::get(quad_vector);
        break;
      case 65:
        internal::LebedevData<65>::get(quad_vector);
        break;
      case 71:
        internal::LebedevData<71>::get(quad_vector);
        break;
      case 77:
        internal::LebedevData<77>::get(quad_vector);
        break;
      case 83:
        internal::LebedevData<83>::get(quad_vector);
        break;
      case 89:
        internal::LebedevData<89>::get(quad_vector);
        break;
      case 95:
        internal::LebedevData<95>::get(quad_vector);
        break;
      case 101:
        internal::LebedevData<101>::get(quad_vector);
        break;
      case 107:
        internal::LebedevData<107>::get(quad_vector);
        break;
      case 113:
        internal::LebedevData<113>::get(quad_vector);
        break;
      case 119:
        internal::LebedevData<119>::get(quad_vector);
        break;
      case 125:
        internal::LebedevData<125>::get(quad_vector);
        break;
      case 131:
        internal::LebedevData<131>::get(quad_vector);
        break;
      default:
        DUNE_THROW(NotImplemented, "Requested order is not available!");
    }

    Dune::QuadratureRule<FieldType, 3> quad_rule;
    for (const auto& quad_point : quad_vector) {
      const auto& pos = quad_point.first;
      quad_rule.emplace_back(FieldVector<FieldType, 3>{{pos[0], pos[1], pos[2]}}, quad_point.second);
    }

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
