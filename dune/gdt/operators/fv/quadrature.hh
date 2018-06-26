// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_QUADRATURE_HH
#define DUNE_GDT_OPERATORS_FV_QUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class DomainFieldType, size_t dimDomain>
struct product_quadrature_helper;

template <class DomainFieldType>
struct product_quadrature_helper<DomainFieldType, 1>
{
  static QuadratureRule<DomainFieldType, 0> get(const QuadratureRule<DomainFieldType, 1>& /*quadrature_1d*/)
  {
    QuadratureRule<DomainFieldType, 0> ret;
    ret.push_back(QuadraturePoint<DomainFieldType, 0>({}, 1));
    return ret;
  }
};

template <class DomainFieldType>
struct product_quadrature_helper<DomainFieldType, 2>
{
  static QuadratureRule<DomainFieldType, 1> get(const QuadratureRule<DomainFieldType, 1>& quadrature_1d)
  {
    return quadrature_1d;
  }
};

template <class DomainFieldType>
struct product_quadrature_helper<DomainFieldType, 3>
{
  static QuadratureRule<DomainFieldType, 2> get(const QuadratureRule<DomainFieldType, 1>& quadrature_1d)
  {
    QuadratureRule<DomainFieldType, 2> ret;
    for (size_t ii = 0; ii < quadrature_1d.size(); ++ii)
      for (size_t jj = 0; jj < quadrature_1d.size(); ++jj)
        ret.push_back(
            QuadraturePoint<DomainFieldType, 2>({quadrature_1d[ii].position()[0], quadrature_1d[jj].position()[0]},
                                                quadrature_1d[ii].weight() * quadrature_1d[jj].weight()));
    return ret;
  }
};


} // namespace internal


template <class DomainFieldType>
QuadratureRule<DomainFieldType, 1> default_1d_quadrature(const size_t reconstructionOrder)
{
  assert(reconstructionOrder <= std::numeric_limits<int>::max());
  return QuadratureRules<DomainFieldType, 1>::rule(Dune::GeometryType(Dune::GeometryType::BasicType::cube, 1),
                                                   static_cast<int>(reconstructionOrder));
}

template <class DomainFieldType, size_t dimDomain>
Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
product_quadrature_on_intersection(const Dune::QuadratureRule<DomainFieldType, 1>& quadrature_1d)
{
  return internal::product_quadrature_helper<DomainFieldType, dimDomain>::get(quadrature_1d);
}

template <class DomainFieldType, size_t dimDomain>
const QuadratureRule<DomainFieldType, dimDomain - 1>& midpoint_quadrature()
{
  static auto midpoint_quadrature_ =
      product_quadrature_on_intersection<DomainFieldType, dimDomain>(default_1d_quadrature<DomainFieldType>(0));
  return midpoint_quadrature_;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_QUADRATURE_HH
