// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/orthonormal.hh>

#include "interfaces.hh"
#include "default.hh"
#include "wrapper.hh"
#include "0d.hh"

namespace Dune {
namespace GDT {


/**
 * \note Update this class if anything changes in dune-localfunctions.
 * \todo Add dimension-and-geometry-type-dependent checks for the order.
 */
template <class D, size_t d, class R>
class LocalOrthonormalFiniteElementFactory
{
  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, 1>;

  template <size_t d_ = d, bool anything = true>
  struct helper
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      return std::make_unique<LocalFiniteElementWrapper<OrthonormalLocalFiniteElement<d, D, R>, D, d, R, 1>>(
          order, geometry_type, order);
    }
  };

  template <bool anything>
  struct helper<0, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& /*order*/)
    {
      // If we need this, and geometry_type.dim() == 0, we must simply implement the corresponding ctors of the 0d FE!
      DUNE_THROW_IF(
          geometry_type.dim() != 0 || !geometry_type.isSimplex(),
          Exceptions::finite_element_error,
          "when creating a local 0d orthonomal finite element: not available for geometry_type = " << geometry_type);
      return std::make_unique<Local0dFiniteElement<D, R>>();
    }
  };

public:
  static std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>> create(const GeometryType& geometry_type,
                                                                         const int& order)
  {
    return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>>(helper<>::create(geometry_type, order));
  }
}; // class LocalOrthonormalFiniteElementFactory


template <class D, size_t d, class R>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>>
make_local_orthonormal_finite_element(const GeometryType& geometry_type, const int order)
{
  return LocalOrthonormalFiniteElementFactory<D, d, R>::create(geometry_type, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_ORTHONORMAL_HH
