// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

/**
 * There is similar functionality in dune/localfunctions/meta/power.hh, but that one is only implemented for the
 * "global" finite elements, while we use the local finite elements.
 **/


#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_POWER_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_POWER_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace GDT {


/**
 * \brief The basis for LocalPowerFiniteElement.
 *
 * \note Only implemented for scalar and vector valued local finite elements.
 *
 * \sa LocalFiniteElementBasisInterface
 * \sa LocalPowerFiniteElement
 */
template <size_t power, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalPowerFiniteElementBasis
{
  static_assert(AlwaysFalse<D>::value, "Not implemented for matrix-valued finite elements");
};


template <size_t power, class D, size_t d, class R, size_t r>
class LocalPowerFiniteElementBasis<power, D, d, R, r, 1> : public LocalFiniteElementBasisInterface<D, d, R, power * r>
{
  using ThisType = LocalPowerFiniteElementBasis;
  using BaseType = LocalFiniteElementBasisInterface<D, d, R, power * r>;

public:
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  using UnpoweredType = LocalFiniteElementBasisInterface<D, d, R, r>;

  LocalPowerFiniteElementBasis(const UnpoweredType& unpowered)
    : unpowered_(unpowered.copy())
  {}

  LocalPowerFiniteElementBasis(const ThisType& other)
    : unpowered_(other.unpowered_->copy())
  {}

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return unpowered_->geometry_type();
  }

  int order() const override final
  {
    return unpowered_->order();
  }

  size_t size() const override final
  {
    return unpowered_->size() * power;
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& point_in_reference_element, std::vector<RangeType>& result) const override final
  {
    const size_t unpowered_sz = unpowered_->size();
    unpowered_->evaluate(point_in_reference_element, unpowered_values_);
    assert(unpowered_values_.size() >= unpowered_sz);
    const size_t sz = this->size();
    if (result.size() < sz)
      result.resize(sz);
    for (size_t pp = 0; pp < power; ++pp)
      for (size_t ii = 0; ii < unpowered_sz; ++ii) {
        result[pp * unpowered_sz + ii] *= 0.;
        for (size_t rr = 0; rr < r; ++rr)
          result[pp * unpowered_sz + ii][pp * r + rr] = unpowered_values_[ii][rr];
      }
  } // ... evaluate(...)

  using BaseType::jacobian;

  void jacobian(const DomainType& point_in_reference_element,
                std::vector<DerivativeRangeType>& result) const override final
  {
    const size_t unpowered_sz = unpowered_->size();
    unpowered_->jacobian(point_in_reference_element, unpowered_jacobians_);
    assert(unpowered_jacobians_.size() >= unpowered_sz);
    const size_t sz = this->size();
    if (result.size() < sz)
      result.resize(sz);
    for (size_t pp = 0; pp < power; ++pp)
      for (size_t ii = 0; ii < unpowered_sz; ++ii) {
        result[pp * unpowered_sz + ii] *= 0.;
        for (size_t rr = 0; rr < r; ++rr)
          result[pp * unpowered_sz + ii][pp * r + rr] = unpowered_jacobians_[ii][rr];
      }
  } // ... jacobian(...)

private:
  const std::unique_ptr<const UnpoweredType> unpowered_;
  mutable std::vector<typename UnpoweredType::RangeType> unpowered_values_;
  mutable std::vector<typename UnpoweredType::DerivativeRangeType> unpowered_jacobians_;
}; // class LocalPowerFiniteElementBasis


/**
 * \brief The interpolation for LocalPowerFiniteElement.
 *
 * \note Only implemented for scalar and vector valued local finite elements.
 *
 * \sa LocalFiniteElementInterpolationInterface
 * \sa LocalPowerFiniteElement
 */
template <size_t power, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalPowerFiniteElementInterpolation
{
  static_assert(AlwaysFalse<D>::value, "Not implemented for matrix-valued finite elements");
};


template <size_t power, class D, size_t d, class R, size_t r>
class LocalPowerFiniteElementInterpolation<power, D, d, R, r, 1>
  : public LocalFiniteElementInterpolationInterface<D, d, R, power * r, 1>
{
  using ThisType = LocalPowerFiniteElementInterpolation;
  using BaseType = LocalFiniteElementInterpolationInterface<D, d, R, power * r, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  using UnpoweredType = LocalFiniteElementInterpolationInterface<D, d, R, r, 1>;

  LocalPowerFiniteElementInterpolation(const UnpoweredType& unpowered)
    : unpowered_(unpowered.copy())
  {}

  LocalPowerFiniteElementInterpolation(const ThisType& other)
    : unpowered_(other.unpowered_->copy())
  {}

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return unpowered_->geometry_type();
  }

  size_t size() const override final
  {
    return unpowered_->size() * power;
  }

  using BaseType::interpolate;

  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   const int order,
                   DynamicVector<R>& dofs) const override final
  {
    const size_t unpowered_sz = unpowered_->size();
    if (unpowered_dofs_.size() < unpowered_sz)
      unpowered_dofs_.resize(unpowered_sz);
    const size_t sz = this->size();
    if (dofs.size() < sz)
      dofs.resize(sz);
    eval_cache_.clear();
    for (size_t pp = 0; pp < power; ++pp) {
      unpowered_->interpolate(
          [&](const auto& point_in_reference_element) {
            auto cache_it = eval_cache_.find(point_in_reference_element);
            if (cache_it == eval_cache_.end())
              cache_it =
                  eval_cache_
                      .insert(std::make_pair(point_in_reference_element, local_function(point_in_reference_element)))
                      .first;
            const RangeType& tmp = cache_it->second;
            FieldVector<R, r> ret;
            for (size_t rr = 0; rr < r; ++rr)
              ret[rr] = tmp[pp * r + rr];
            return ret;
          },
          order,
          unpowered_dofs_);
      assert(unpowered_dofs_.size() >= unpowered_sz);
      for (size_t ii = 0; ii < unpowered_sz; ++ii)
        dofs[pp * unpowered_sz + ii] = unpowered_dofs_[ii];
    }
  } // ... interpolate(...)

private:
  const std::unique_ptr<const UnpoweredType> unpowered_;
  mutable DynamicVector<R> unpowered_dofs_;
  mutable std::map<DomainType, RangeType, XT::Common::FieldVectorLess> eval_cache_;
}; // class LocalPowerFiniteElementInterpolation


/**
 * \brief The coefficients for LocalPowerFiniteElement.
 *
 * \note Only implemented for scalar and vector valued local finite elements.
 *
 * \sa LocalFiniteElementCoefficientsInterface
 * \sa LocalPowerFiniteElement
 */
template <class D, size_t d>
class LocalPowerFiniteElementCoefficients : public LocalFiniteElementCoefficientsInterface<D, d>
{
  using ThisType = LocalPowerFiniteElementCoefficients;
  using BaseType = LocalFiniteElementCoefficientsInterface<D, d>;

public:
  using UnpoweredType = BaseType;

  LocalPowerFiniteElementCoefficients(const UnpoweredType& unpowered, const size_t power)
    : unpowered_(unpowered.copy())
    , power_(power)
    , local_keys_(unpowered_->size() * power_)
  {
    const size_t unpowered_sz = unpowered_->size();
    const auto unpowered_indices = unpowered_->local_key_indices();
    for (size_t pp = 0; pp < power; ++pp)
      for (size_t ii = 0; ii < unpowered_sz; ++ii) {
        const auto& unpowered_local_key = unpowered_->local_key(ii);
        const auto sub_entity = unpowered_local_key.subEntity();
        const auto codim = unpowered_local_key.codim();
        const auto num_unpowered_indices = unpowered_indices[codim][sub_entity].size();
        local_keys_[pp * unpowered_sz + ii] =
            LocalKey(sub_entity,
                     codim,
                     XT::Common::numeric_cast<unsigned int>(pp * num_unpowered_indices + unpowered_local_key.index()));
      }
  } // LocalPowerFiniteElementCoefficients(...)

  LocalPowerFiniteElementCoefficients(const ThisType& other)
    : unpowered_(other.unpowered_->copy())
    , power_(other.power_)
    , local_keys_(other.local_keys_)
  {}

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return unpowered_->geometry_type();
  }

  size_t size() const override final
  {
    return unpowered_->size() * power_;
  }

  const LocalKey& local_key(const size_t ii) const override final
  {
    if (ii >= size())
      DUNE_THROW(Exceptions::finite_element_error, "ii = " << ii << "\n   size() = " << size());
    return local_keys_[ii];
  }

private:
  const std::unique_ptr<const UnpoweredType> unpowered_;
  const size_t power_;
  std::vector<LocalKey> local_keys_;
}; // class LocalPowerFiniteElementCoefficients


/**
 * \brief Models the product of a given local finite element with itself, power times.
 *
 * \note Only implemented for scalar and vector valued local finite elements.
 *
 * \sa LocalFiniteElementInterface
 * \sa LocalPowerFiniteElementBasis
 * \sa LocalPowerFiniteElementInterpolation
 * \sa LocalPowerFiniteElementCoefficients
 */
template <size_t power, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalPowerFiniteElement
{
  static_assert(AlwaysFalse<D>::value, "Not implemented for matrix-valued finite elements");
};


// Specialization for the unpowered case.
template <class D, size_t d, class R, size_t r, size_t rC>
class LocalPowerFiniteElement<1, D, d, R, r, rC> : public LocalFiniteElementDefault<D, d, R, r, rC>
{
  using ThisType = LocalPowerFiniteElement;
  using BaseType = LocalFiniteElementDefault<D, d, R, r, rC>;

public:
  using typename BaseType::DomainType;

  using UnpoweredType = LocalFiniteElementInterface<D, d, R, r, rC>;

  LocalPowerFiniteElement(const UnpoweredType& unpowered)
    : BaseType(unpowered.order(),
               unpowered.basis().copy(),
               unpowered.coefficients().copy(),
               unpowered.interpolation().copy(),
               unpowered.is_lagrangian() ? unpowered.lagrange_points() : std::vector<DomainType>())
  {}
}; // class LocalPowerFiniteElement<1, ...>


// The actual implementation, but only for scalar and vector valued FEs.
template <size_t power, class D, size_t d, class R, size_t r>
class LocalPowerFiniteElement<power, D, d, R, r, 1> : public LocalFiniteElementDefault<D, d, R, power * r>
{
  using ThisType = LocalPowerFiniteElement;
  using BaseType = LocalFiniteElementDefault<D, d, R, power * r>;

public:
  using typename BaseType::DomainType;

  using UnpoweredType = LocalFiniteElementInterface<D, d, R, r>;

  LocalPowerFiniteElement(const UnpoweredType& unpowered)
    : BaseType(unpowered.order(),
               new LocalPowerFiniteElementBasis<power, D, d, R, r>(unpowered.basis()),
               new LocalPowerFiniteElementCoefficients<D, d>(unpowered.coefficients(), power),
               new LocalPowerFiniteElementInterpolation<power, D, d, R, r>(unpowered.interpolation()),
               unpowered.is_lagrangian() ? unpowered.lagrange_points() : std::vector<DomainType>())
  {}
}; // class LocalPowerFiniteElement


template <size_t power, class D, size_t d, class R, size_t r>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, power * r>>
make_local_powered_finite_element(const LocalFiniteElementInterface<D, d, R, r>& unpowered_finite_element)
{
  return std::make_unique<LocalPowerFiniteElement<power, D, d, R, r>>(unpowered_finite_element);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_POWER_HH
