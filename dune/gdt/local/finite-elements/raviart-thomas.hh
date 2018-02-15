// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH

#include <dune/localfunctions/raviartthomas.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"
#include "wrapper.hh"

namespace Dune {
namespace GDT {
namespace internal {


/**
 * \note Update this class if anything changes in dune-localfunctions.
 */
template <class D, size_t d, class R>
class RaviartThomasLocalFiniteElementFactory
{
  static_assert(1 <= d && d <= 3, "There is no local Raviart-Thomas finite element available for other dimension!");

  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, d, 1>;

  template <size_t d_ = d, bool anything = true>
  struct helper;

  template <bool anything>
  struct helper<1, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& polorder)
    {
      // everything is a simplex in 1d
      using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
      return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE(geometry_type, polorder));
    }
  };

  template <bool anything>
  struct helper<2, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& polorder)
    {
      if (geometry_type.isSimplex()) {
        using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
        return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE(geometry_type, polorder));
      } else if (geometry_type.isCube()) {
        if (polorder == 0) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 0>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else if (polorder == 1) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 1>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else if (polorder == 2) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 2>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else if (polorder == 3) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 3>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else if (polorder == 4) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 4>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else
          DUNE_THROW(
              Exceptions::finite_element_error,
              "when creating a local Raviart-Thomas finite element: there is none available on cubes for polorder "
                  << polorder
                  << " (if you think there is, update this class)!");
      } else
        DUNE_THROW(Exceptions::finite_element_error,
                   "when creating a local Raviart-Thomas finite element: there is none available for the "
                   "following geometry type: "
                       << geometry_type
                       << "(if you think there is, update this class)!");
    }
  }; // helper<2, ...>

  template <bool anything>
  struct helper<3, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& polorder)
    {
      if (geometry_type.isSimplex()) {
        using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
        return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE(geometry_type, polorder));
      } else if (geometry_type.isCube()) {
        if (polorder == 0) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 0>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else if (polorder == 1) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 1>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(new FE());
        } else
          DUNE_THROW(
              Exceptions::finite_element_error,
              "when creating a local Raviart-Thomas finite element: there is none available on cubes for polorder "
                  << polorder
                  << " (if you think there is, update this class)!");
      } else
        DUNE_THROW(Exceptions::finite_element_error,
                   "when creating a local Raviart-Thomas finite element: there is none available for the "
                   "following geometry type: "
                       << geometry_type
                       << "(if you think there is, update this class)!");
    }
  }; // helper<3, ...>

public:
  static std::unique_ptr<LocalFiniteElementInterface<D, d, R, d, 1>> create(const GeometryType& geometry_type,
                                                                            const int& polorder)
  {
    return std::unique_ptr<LocalFiniteElementInterface<D, d, R, d, 1>>(helper<>::create(geometry_type, polorder));
  }
}; // class RaviartThomasLocalFiniteElementFactory


} // namespace internal


template <class D, size_t d, class R>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, d, 1>>
make_raviart_thomas_local_finite_element(const GeometryType& geometry_type, const int& polorder)
{
  return internal::RaviartThomasLocalFiniteElementFactory<D, d, R>::create(geometry_type, polorder);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH
