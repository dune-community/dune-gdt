// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_LAGRANGE_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_LAGRANGE_HH

#include <memory>
#include <sstream>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"
#include "default.hh"
#include "power.hh"
#include "wrapper.hh"
#include "0d.hh"

namespace Dune {
namespace GDT {


/**
 * Wraps the P0LocalFiniteElement from dune-localfunctions and adds Lagrange points.
 */
template <class D, size_t d, class R>
class LocalZeroOrderLagrangeFiniteElement
    : XT::Common::ConstStorageProvider<LocalFiniteElementWrapper<P0LocalFiniteElement<D, R, d>, D, d, R, 1>>,
      public LocalFiniteElementDefault<D, d, R, 1>
{
  using ThisType = LocalZeroOrderLagrangeFiniteElement<D, d, R>;
  using Implementation = P0LocalFiniteElement<D, R, d>;
  using Wrapper = LocalFiniteElementWrapper<Implementation, D, d, R, 1>;
  using Storage = XT::Common::ConstStorageProvider<Wrapper>;
  using BaseType = LocalFiniteElementDefault<D, d, R, 1>;

public:
  LocalZeroOrderLagrangeFiniteElement(const GeometryType& geometry_type)
    : XT::Common::ConstStorageProvider<Wrapper>(new Wrapper(0, geometry_type))
    , BaseType(0,
               Storage::access().basis().copy(),
               Storage::access().coefficients().copy(),
               Storage::access().interpolation().copy(),
               {ReferenceElements<D, d>::general(geometry_type).position(0, 0)})
  {
  }
}; // class LocalZeroOrderLagrangeFiniteElement


/**
 * \note Update this class if anything changes in dune-localfunctions.
 */
template <class D, size_t d, class R, size_t r = 1>
class LocalLagrangeFiniteElementFactory
{
  using ScalarLocalFiniteElementType = LocalFiniteElementInterface<D, d, R, 1>;
  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, r>;

  static std::string order_error(const GeometryType& geometry_type, const int order)
  {
    std::stringstream ss;
    ss << "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in " << d
       << "d on a " << geometry_type << " reference element for order " << order
       << " (if you think it is working, update this check)!";
    return ss.str();
  }

  static std::string geometry_error(const GeometryType& geometry_type, const int order)
  {
    std::stringstream ss;
    ss << "when creating a local Lagrange finite element: this is untested!\n"
       << "Please update this check if you believe that LagrangeLocalFiniteElement is available for\n- dimension: " << d
       << "\n- geometry_type: " << geometry_type << "\n- order: " << order;
    return ss.str();
  }

  // Fist we create the scalar FE ...

  template <size_t d_ = d, bool anything = true>
  struct scalar_helper
  {
    static std::unique_ptr<ScalarLocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      DUNE_THROW_IF(geometry_type.dim() != d,
                    Exceptions::finite_element_error,
                    "geometry_type = " << geometry_type << "\nd = " << d);
      // special case
      if (order == 0)
        return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>>(
            new LocalZeroOrderLagrangeFiniteElement<D, d, R>(geometry_type));
      // checks
      if (d == 1)
        DUNE_THROW_IF(order > 18, Exceptions::finite_element_error, order_error(geometry_type, order));
      else if (d == 2) {
        if (geometry_type == GeometryType(GeometryType::simplex, 2))
          DUNE_THROW_IF(order > 15, Exceptions::finite_element_error, order_error(geometry_type, order));
        else if (geometry_type == GeometryType(GeometryType::cube, 2))
          DUNE_THROW_IF(order > 10, Exceptions::finite_element_error, order_error(geometry_type, order));
        else
          DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      } else if (d == 3) {
        if (geometry_type == GeometryType(GeometryType::simplex, 3))
          DUNE_THROW_IF(order > 14, Exceptions::finite_element_error, order_error(geometry_type, order));
        else if (geometry_type == GeometryType(GeometryType::cube, 3))
          DUNE_THROW_IF(order > 7, Exceptions::finite_element_error, order_error(geometry_type, order));
        else if (geometry_type == GeometryType(GeometryType::prism, 3))
          DUNE_THROW_IF(order > 9, Exceptions::finite_element_error, order_error(geometry_type, order));
        else if (geometry_type == GeometryType(GeometryType::pyramid, 3))
          DUNE_THROW(
              Exceptions::finite_element_error,
              "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail on "
              "pyramids (if you think it is working, update this check)!");
        else
          DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      } else if (d == 4) {
        if (geometry_type == GeometryType(GeometryType::simplex, 4))
          DUNE_THROW_IF(order > 13, Exceptions::finite_element_error, order_error(geometry_type, order));
        else if (geometry_type == GeometryType(GeometryType::cube, 4))
          DUNE_THROW_IF(order > 5, Exceptions::finite_element_error, order_error(geometry_type, order));
        else
          DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      } else
        DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      // the actual finite element
      return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>>(
          new LocalFiniteElementWrapper<LagrangeLocalFiniteElement<EquidistantPointSet, d, D, R>, D, d, R, 1>(
              order, geometry_type, order));
    }
  }; // helper<...>

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
  }; // helper<...>

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
}; // class LocalLagrangeFiniteElementFactory


template <class D, size_t d, class R, size_t r = 1>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, r>>
make_local_lagrange_finite_element(const GeometryType& geometry_type, const int order)
{
  return LocalLagrangeFiniteElementFactory<D, d, R, r>::create(geometry_type, order);
}


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_LAGRANGE_HH
