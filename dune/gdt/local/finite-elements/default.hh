// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH

#include <dune/xt/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \brief Implements LocalFiniteElementInterface, given a basis, coefficients and an interpolation.
 */
template <class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementDefault : public LocalFiniteElementInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementDefault<D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterface<D, d, R, r, rC>;

public:
  using typename BaseType::BasisType;
  using typename BaseType::CoefficientsType;
  using typename BaseType::DomainType;
  using typename BaseType::InterpolationType;

  LocalFiniteElementDefault(const int ord,
                            const BasisType& bas,
                            const CoefficientsType& coeffs,
                            const InterpolationType& inter,
                            std::vector<DomainType> lps = {})
    : geometry_type_(bas.geometry_type())
    , order_(ord)
    , basis_(bas)
    , coefficients_(coeffs)
    , interpolation_(inter)
    , lagrange_points_(lps)
  {
    check_input();
  }

  /**
   * \attention Do not delete any of the moved in raw pointers manually afterwads!
   **/
  LocalFiniteElementDefault(const int ord,
                            BasisType*&& bas_ptr,
                            CoefficientsType*&& coeffs_ptr,
                            InterpolationType*&& inter_ptr,
                            std::vector<DomainType> lps = {})
    : geometry_type_(bas_ptr->geometry_type())
    , order_(ord)
    , basis_(std::move(bas_ptr))
    , coefficients_(std::move(coeffs_ptr))
    , interpolation_(std::move(inter_ptr))
    , lagrange_points_(lps)
  {
    check_input();
  }

  const GeometryType& geometry_type() const override
  {
    return geometry_type_;
  }

  int order() const override
  {
    return order_;
  }

  size_t size() const override
  {
    return basis_.access().size();
  }

  const BasisType& basis() const override
  {
    return basis_.access();
  }

  const CoefficientsType& coefficients() const override
  {
    return coefficients_.access();
  }

  const InterpolationType& interpolation() const override
  {
    return interpolation_.access();
  }

  bool is_lagrangian() const override
  {
    return lagrange_points_.size() > 0;
  }

  const std::vector<DomainType>& lagrange_points() const override
  {
    if (!is_lagrangian())
      DUNE_THROW(Exceptions::finite_element_error, "do not call lagrange_points() if is_lagrangian() is false!");
    return lagrange_points_;
  }

private:
  void check_input()
  {
    DUNE_THROW_IF(coefficients_.access().geometry_type() != geometry_type_,
                  Exceptions::finite_element_error,
                  "\n   coefficients_.access().geometry_type() = " << coefficients_.access().geometry_type()
                                                                   << "\n   geometry_type_ = " << geometry_type_);
    DUNE_THROW_IF(coefficients_.access().size() != basis_.access().size(),
                  Exceptions::finite_element_error,
                  "\n   coefficients_.access().size() = "
                      << coefficients_.access().size() << "\n   basis_.access().size() = " << basis_.access().size());
    DUNE_THROW_IF(interpolation_.access().geometry_type() != geometry_type_,
                  Exceptions::finite_element_error,
                  "\n   interpolation_.access().geometry_type() = " << interpolation_.access().geometry_type()
                                                                    << "\n   " << geometry_type_);
    DUNE_THROW_IF(interpolation_.access().size() != basis_.access().size(),
                  Exceptions::finite_element_error,
                  "\n   interpolation_.access().size() = "
                      << interpolation_.access().size() << "\n   basis_.access().size() = " << basis_.access().size());
    DUNE_THROW_IF(!(lagrange_points_.size() == basis_.access().size() / (r * rC) || lagrange_points_.size() == 0),
                  Exceptions::finite_element_error,
                  "\n   lagrange_points_.size() = " << lagrange_points_.size()
                                                    << "\n   basis_.access().size() / (r * rC) = "
                                                    << basis_.access().size() / (r * rC));
  } // ... check_input(...)

  const GeometryType geometry_type_;
  const int order_;
  const XT::Common::ConstStorageProvider<BasisType> basis_;
  const XT::Common::ConstStorageProvider<CoefficientsType> coefficients_;
  const XT::Common::ConstStorageProvider<InterpolationType> interpolation_;
  const std::vector<DomainType> lagrange_points_;
}; // class LocalFiniteElementDefault


template <class D, size_t d, class R = double, size_t r = 1, size_t rC = 1>
class ThreadSafeDefaultLocalFiniteElementFamily : public LocalFiniteElementFamilyInterface<D, d, R, r, rC>
{
  using ThisType = ThreadSafeDefaultLocalFiniteElementFamily<D, d, R, r, rC>;
  using BaseType = LocalFiniteElementFamilyInterface<D, d, R, r, rC>;

public:
  using typename BaseType::LocalFiniteElementType;

  ThreadSafeDefaultLocalFiniteElementFamily(
      std::function<std::unique_ptr<LocalFiniteElementType>(const GeometryType&, const int&)> factory)
    : factory_(factory)
  {}

  ThreadSafeDefaultLocalFiniteElementFamily(const ThisType& other)
    : factory_(other.factory_)
  {
    // we do not even try to create the FEs in a thread safe way, they will just be recreated when required
  }

  ThreadSafeDefaultLocalFiniteElementFamily(ThisType&& source)
    : factory_(std::move(source.factory_))
    , fes_(std::move(source.fes_))
  {}

  const LocalFiniteElementType& get(const GeometryType& geometry_type, const int order) const override final
  {
    const auto key = std::make_pair(geometry_type, order);
    // if the FE already exists, no need to lock since at() is thread safe and we are returning the object reference,
    // not a map iterator which might get invalidated
    // TODO: Double checked locking pattern is not thread-safe without memory barriers.
    if (fes_.count(key) == 0) {
      // the FE needs to be created, we need to lock
      std::lock_guard<std::mutex> DXTC_UNUSED(guard)(mutex_);
      // and to check again if someone else created the FE while we were waiting to acquire the lock
      if (fes_.count(key) == 0) {
        auto dings = factory_(geometry_type, order);
        fes_.insert(std::make_pair(key, std::shared_ptr<LocalFiniteElementType>(dings.release())));
      }
    }
    return *fes_.at(key);
  } // ... get(...)

private:
  const std::function<std::unique_ptr<LocalFiniteElementType>(const GeometryType&, const int&)> factory_;
  mutable std::map<std::pair<GeometryType, int>, std::shared_ptr<LocalFiniteElementType>> fes_;
  mutable std::mutex mutex_;
}; // class ThreadSafeDefaultLocalFiniteElementFamily


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH
