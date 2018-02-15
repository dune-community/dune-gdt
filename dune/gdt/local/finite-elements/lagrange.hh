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

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"
#include "wrapper.hh"

namespace Dune {
namespace GDT {


/**
 * This is the P0LocalFiniteElement from dune-localfunctions with suitable Lagrange points (which is also all that
 * differs in the implementation with respect to LocalFiniteElementWrapper).
 */
template <class D, size_t d, class R>
class P0LagrangeFiniteElement : public LocalFiniteElementInterface<D, d, R, 1, 1>
{
  using ThisType = P0LagrangeFiniteElement<D, d, R>;
  using BaseType = LocalFiniteElementInterface<D, d, R, 1, 1>;

  using Implementation = P0LocalFiniteElement<D, R, d>;
  using BasisWrapperType =
      LocalFiniteElementBasisWrapper<typename Implementation::Traits::LocalBasisType, D, d, R, 1, 1>;
  using CoefficientsWrapperType =
      LocalFiniteElementCoefficientsWrapper<typename Implementation::Traits::LocalCoefficientsType>;
  using InterpolationWrapperType =
      LocalFiniteElementInterpolationWrapper<typename Implementation::Traits::LocalInterpolationType, D, d, R, 1, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::BasisType;
  using typename BaseType::CoefficientsType;
  using typename BaseType::InterpolationType;

  P0LagrangeFiniteElement(Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_({ReferenceElements<D, d>::general(imp_.access().type()).position(0, 0)})
  {
  }

  P0LagrangeFiniteElement(const Implementation& imp)
    : imp_(imp)
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_({ReferenceElements<D, d>::general(imp_.access().type()).position(0, 0)})
  {
  }

  template <class... Args>
  explicit P0LagrangeFiniteElement(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_({ReferenceElements<D, d>::general(imp_.access().type()).position(0, 0)})
  {
  }

  GeometryType geometry_type() const
  {
    return imp_.access().type();
  }

  size_t size() const override final
  {
    return XT::Common::numeric_cast<size_t>(imp_.access().size());
  }

  const BasisType& basis() const override final
  {
    return basis_;
  }

  const CoefficientsType& coefficients() const override final
  {
    return coefficients_;
  }

  const InterpolationType& interpolation() const override final
  {
    return interpolation_;
  }

  bool is_lagrangian() const override final
  {
    return true;
  }

  const std::vector<DomainType>& lagrange_points() const override final
  {
    return lagrange_points_;
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
  const BasisWrapperType basis_;
  const CoefficientsWrapperType coefficients_;
  const InterpolationWrapperType interpolation_;
  const std::vector<DomainType> lagrange_points_;
}; // class P0LagrangeFiniteElement


template <class D, size_t d, class R>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1, 1>>
make_lagrange_local_finite_element(const GeometryType& geometry_type, const int& polorder)
{
  // special case
  if (d > 0 && polorder == 0)
    return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1, 1>>(
        new P0LagrangeFiniteElement<D, d, R>(geometry_type));
  // checks
  if (d == 1) {
    if (polorder > 18)
      DUNE_THROW(Exceptions::finite_element_error,
                 "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 1d "
                 "for polorder "
                     << polorder
                     << " (if you think it is working, update this check)!");
  } else if (d == 2) {
    if (geometry_type == GeometryType(GeometryType::simplex, 2)) {
      if (polorder > 15)
        DUNE_THROW(
            Exceptions::finite_element_error,
            "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 2d "
            "on simplices for polorder "
                << polorder
                << " (if you think it is working, update this check)!");
    } else if (geometry_type == GeometryType(GeometryType::cube, 2)) {
      if (polorder > 10)
        DUNE_THROW(
            Exceptions::finite_element_error,
            "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 2d "
            "on cubes for polorder "
                << polorder
                << " (if you think it is working, update this check)!");
    } else
      DUNE_THROW(Exceptions::finite_element_error,
                 "when creating a local Lagrange finite element: this is untested!\n"
                     << "Please update this check if you believe that a suitable finite element is available for\n"
                     << "- dimension: "
                     << d
                     << "\n"
                     << "-n geometry_type: "
                     << geometry_type
                     << "\n"
                     << "- polorder: "
                     << polorder);
  } else if (d == 3) {
    if (geometry_type == GeometryType(GeometryType::simplex, 3)) {
      if (polorder > 14)
        DUNE_THROW(
            Exceptions::finite_element_error,
            "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 3d "
            "on simplices for polorder "
                << polorder
                << " (if you think it is working, update this check)!");
    } else if (geometry_type == GeometryType(GeometryType::cube, 3)) {
      if (polorder > 7)
        DUNE_THROW(
            Exceptions::finite_element_error,
            "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 3d "
            "on cubes for polorder "
                << polorder
                << " (if you think it is working, update this check)!");
    } else if (geometry_type == GeometryType(GeometryType::prism, 3)) {
      if (polorder > 9)
        DUNE_THROW(
            Exceptions::finite_element_error,
            "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 3d "
            "on prisms for polorder "
                << polorder
                << " (if you think it is working, update this check)!");
    } else if (geometry_type == GeometryType(GeometryType::pyramid, 3)) {
      DUNE_THROW(Exceptions::finite_element_error,
                 "when creating a local Lagrange finite element: the LagrangeLocalFiniteElement is known to fail in 3d "
                 "on pyramids for polorder "
                     << polorder
                     << " (if you think it is working, update this check)!");
    } else
      DUNE_THROW(Exceptions::finite_element_error,
                 "when creating a local Lagrange finite element: this is untested!\n"
                     << "Please update this check if you believe that a suitable finite element is available for\n"
                     << "- dimension: "
                     << d
                     << "\n"
                     << "-n geometry_type: "
                     << geometry_type
                     << "\n"
                     << "- polorder: "
                     << polorder);
  } else
    DUNE_THROW(Exceptions::finite_element_error,
               "when creating a local Lagrange finite element: this is untested!\n"
                   << "Please update this check if you believe that a suitable finite element is available for\n"
                   << "- dimension: "
                   << d
                   << "\n"
                   << "-n geometry_type: "
                   << geometry_type
                   << "\n"
                   << "- polorder: "
                   << polorder);
  // the actual finite element
  return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1, 1>>(
      new LocalFiniteElementWrapper<LagrangeLocalFiniteElement<EquidistantPointSet, d, D, R>, D, d, R, 1>(geometry_type,
                                                                                                          polorder));
} // ... make_lagrange_local_finite_element(...)


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_LAGRANGE_HH
