// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/orthonormal.hh>

#include "interfaces.hh"
#include "default.hh"
#include "power.hh"
#include "wrapper.hh"
#include "0d.hh"

namespace Dune {
namespace GDT {


/**
 * \note Update this class if anything changes in dune-localfunctions.
 * \todo Add dimension-and-geometry-type-dependent checks for the order.
 */
template <class D, size_t d, class R, size_t r = 1>
class LocalOrthonormalFiniteElementFactory
{
  using ScalarLocalFiniteElementType = LocalFiniteElementInterface<D, d, R, 1>;
  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, r>;

  // Fist we create the scalar FE ...

  template <size_t d_ = d, bool anything = true>
  struct scalar_helper // d != 0
  {
    static std::unique_ptr<ScalarLocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      return std::make_unique<LocalFiniteElementWrapper<OrthonormalLocalFiniteElement<d, D, R>, D, d, R, 1>>(
          order, geometry_type, order);
    }
  };

  template <bool anything>
  struct scalar_helper<0, anything>
  {
    static std::unique_ptr<ScalarLocalFiniteElementType> create(const GeometryType& geometry_type, const int& /*order*/)
    {
      // If we need this, and geometry_type.dim() == 0, we must simply implement the corresponding ctors of the 0d FE!
      DUNE_THROW_IF(
          geometry_type.dim() != 0 || !geometry_type.isSimplex(),
          Exceptions::finite_element_error,
          "when creating a local 0d orthonomal finite element: not available for geometry_type = " << geometry_type);
      return std::make_unique<Local0dFiniteElement<D, R>>();
    }
  };

  // ... then we wrap this in a power FE, if required.

  template <size_t d_ = d, size_t r_ = r>
  struct helper // r != 1
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      return make_local_powered_finite_element<r>(*scalar_helper<>::create(geometry_type, order));
    }
  };

  template <size_t d_>
  struct helper<d_, 1>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      return scalar_helper<>::create(geometry_type, order);
    }
  };

public:
  static std::unique_ptr<LocalFiniteElementInterface<D, d, R, r>> create(const GeometryType& geometry_type,
                                                                         const int& order)
  {
    return helper<>::create(geometry_type, order);
  }
}; // class LocalOrthonormalFiniteElementFactory


template <class D, size_t d, class R, size_t r = 1>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, r>>
make_local_orthonormal_finite_element(const GeometryType& geometry_type, const int order)
{
  return LocalOrthonormalFiniteElementFactory<D, d, R, r>::create(geometry_type, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH
