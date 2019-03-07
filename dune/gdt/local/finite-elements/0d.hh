// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_0D_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_0D_HH

#include <dune/geometry/referenceelements.hh>

#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace GDT {


/**
 * \brief The basis for Local0dFiniteElement.
 *
 * \sa LocalFiniteElementBasisInterface
 * \sa Local0dFiniteElement
 */
template <class D, class R>
class Local0dFiniteElementBasis : public LocalFiniteElementBasisInterface<D, 0, R, 1>
{
  using ThisType = Local0dFiniteElementBasis<D, R>;
  using BaseType = LocalFiniteElementBasisInterface<D, 0, R, 1>;

public:
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  Local0dFiniteElementBasis()
    : geometry_type_(GeometryTypes::simplex(0))
  {}

  Local0dFiniteElementBasis(const ThisType& other) = default;

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  int order() const override final
  {
    return 0;
  }

  size_t size() const override final
  {
    return 1;
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& /*point_in_reference_element*/, std::vector<RangeType>& result) const override final
  {
    if (result.size() < 1)
      result.resize(1);
    result[0] = 1.;
  }

  using BaseType::jacobian;

  void jacobian(const DomainType& /*point_in_reference_element*/,
                std::vector<DerivativeRangeType>& result) const override final
  {
    if (result.size() < 1)
      result.resize(1);
    result[0] = 0.;
  }

private:
  const GeometryType geometry_type_;
}; // class Local0dFiniteElementBasis


/**
 * \brief The interpolation for Local0dFiniteElement.
 *
 * \sa LocalFiniteElementInterpolationInterface
 * \sa Local0dFiniteElement
 */
template <class D, class R>
class Local0dFiniteElementInterpolation : public LocalFiniteElementInterpolationInterface<D, 0, R, 1>
{
  using ThisType = Local0dFiniteElementInterpolation<D, R>;
  using BaseType = LocalFiniteElementInterpolationInterface<D, 0, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  Local0dFiniteElementInterpolation()
    : geometry_type_(GeometryTypes::simplex(0))
  {}

  Local0dFiniteElementInterpolation(const ThisType& other) = default;

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  size_t size() const override final
  {
    return 1;
  }

  using BaseType::interpolate;

  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   const int /*order*/,
                   DynamicVector<R>& dofs) const override final
  {
    if (dofs.size() < 1)
      dofs.resize(1);
    dofs[0] = local_function({});
  }

private:
  const GeometryType geometry_type_;
}; // class Local0dFiniteElementInterpolation


/**
 * \brief The coefficients for Local0dFiniteElement.
 *
 * \sa LocalFiniteElementCoefficientsInterface
 * \sa Local0dFiniteElement
 */
template <class D>
class Local0dFiniteElementCoefficients : public LocalFiniteElementCoefficientsInterface<D, 0>
{
  using ThisType = Local0dFiniteElementCoefficients<D>;
  using BaseType = LocalFiniteElementCoefficientsInterface<D, 0>;

public:
  Local0dFiniteElementCoefficients()
    : geometry_type_(GeometryTypes::simplex(0))
    , local_key_(0, 0, 0)
  {}

  Local0dFiniteElementCoefficients(const ThisType& other) = default;

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  size_t size() const override final
  {
    return 1;
  }

  const LocalKey& local_key(const size_t ii) const override final
  {
    DUNE_THROW_IF(ii > 0, Exceptions::finite_element_error, "ii = " << ii << "\n   size() = 1");
    return local_key_;
  }

private:
  const GeometryType geometry_type_;
  const LocalKey local_key_;
}; // class Local0dFiniteElementCoefficients


/**
 * \brief A local finite element on point-like reference elements.
 *
 * We require this, for instance, when we need a basis associated with an intersection (which, in 1d, is a point).
 *
 * \sa LocalFiniteElementInterface
 * \sa Local0dFiniteElementBasis
 * \sa Local0dFiniteElementInterpolation
 * \sa Local0dFiniteElementCoefficients
 */
template <class D, class R>
class Local0dFiniteElement : public LocalFiniteElementDefault<D, 0, R, 1>
{
  using BaseType = LocalFiniteElementDefault<D, 0, R, 1>;

public:
  Local0dFiniteElement()
    : BaseType(0,
               new Local0dFiniteElementBasis<D, R>(),
               new Local0dFiniteElementCoefficients<D>(),
               new Local0dFiniteElementInterpolation<D, R>(),
               {ReferenceElements<D, 0>::simplex().position(0, 0)})
  {}
}; // class Local0dFiniteElement


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_0D_HH
