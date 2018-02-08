// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH

#include <functional>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>

#include <dune/xt/functions/interfaces/local-functions.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


template <class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementBasisInterface
{
  using ThisType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;

public:
  using DomainType = FieldVector<D, d>;
  using RangeType = typename XT::Functions::RangeTypeSelector<R, r, rC>::type;
  using JacobianRangeType = typename XT::Functions::JacobianRangeTypeSelector<d, R, r, rC>::type;

  LocalFiniteElementBasisInterface() = default;
  LocalFiniteElementBasisInterface(const ThisType&) = default;
  LocalFiniteElementBasisInterface(ThisType&&) = default;

  virtual ~LocalFiniteElementBasisInterface() = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  virtual int order() const = 0;

  virtual size_t size() const = 0;

  virtual std::vector<RangeType> evaluate(const DomainType& /*xx*/) const = 0;

  virtual std::vector<JacobianRangeType> jacobian(const DomainType& /*xx*/) const = 0;
}; // class LocalFiniteElementBasisInterface


template <class D, size_t d, class R, size_t r, size_t rC = 1, class F = R>
class LocalFiniteElementInterpolationInterface
{
  using ThisType = LocalFiniteElementInterpolationInterface<D, d, R, r, rC, F>;

public:
  using DomainType = FieldVector<D, d>;
  using RangeType = typename XT::Functions::RangeTypeSelector<R, r, rC>::type;

  LocalFiniteElementInterpolationInterface() = default;
  LocalFiniteElementInterpolationInterface(const ThisType&) = default;
  LocalFiniteElementInterpolationInterface(ThisType&&) = default;

  virtual ~LocalFiniteElementInterpolationInterface() = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  virtual std::vector<F> interpolate(const std::function<RangeType(DomainType)>& local_function) const = 0;
}; // class LocalFiniteElementInterpolationInterface


class LocalFiniteElementCoefficientsInterface
{
  using ThisType = LocalFiniteElementCoefficientsInterface;

public:
  LocalFiniteElementCoefficientsInterface() = default;
  LocalFiniteElementCoefficientsInterface(const ThisType&) = default;
  LocalFiniteElementCoefficientsInterface(ThisType&&) = default;

  virtual ~LocalFiniteElementCoefficientsInterface() = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  virtual size_t size() const = 0;

  virtual const LocalKey& local_key(const size_t ii) const = 0;
}; // class LocalFiniteElementCoefficientsInterface


template <class D, size_t d, class R, size_t r, size_t rC = 1, class F = R>
class LocalFiniteElementInterface
{
  using ThisType = LocalFiniteElementInterface<D, d, R, r, rC, F>;

public:
  using DomainType = FieldVector<D, d>;
  using BasisType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;
  using CoefficientsType = LocalFiniteElementCoefficientsInterface;
  using InterpolationType = LocalFiniteElementInterpolationInterface<D, d, R, r, rC, F>;

  LocalFiniteElementInterface() = default;
  LocalFiniteElementInterface(const ThisType&) = default;
  LocalFiniteElementInterface(ThisType&&) = default;

  virtual ~LocalFiniteElementInterface() = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  virtual GeometryType geometry_type() const = 0;

  virtual size_t size() const = 0;

  virtual const BasisType& basis() const = 0;

  virtual const CoefficientsType& coefficients() const = 0;

  virtual const InterpolationType& interpolation() const = 0;

  virtual bool is_lagrangian() const
  {
    return false;
  }

  virtual const std::vector<DomainType>& lagrange_points() const
  {
    if (is_lagrangian())
      DUNE_THROW(NotImplemented,
                 "the implementor of this local finite element has to provider lagrange_points(), if "
                 "is_lagrangian() is true!");
    else
      DUNE_THROW(finite_element_error, "do not call lagrange_points() if is_lagrangian() is false!");
  }
}; // class LocalFiniteElementInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
