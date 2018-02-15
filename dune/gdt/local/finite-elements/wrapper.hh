// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH

#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace internal {


template <class Implementation, class DomainType>
class LocalFiniteElementInterpolationLagrangepointsAccessor
{
  DXTC_has_method_initialize_once(lagrangePoints);
  static const constexpr bool has_lagrange_points = DXTC_has_method(lagrangePoints)<Implementation>::value;

  template <bool available = has_lagrange_points, bool anything = true>
  struct lagrangian_helper // available == false
  {
    static std::vector<DomainType> get(const Implementation& /*imp*/)
    {
      return {};
    }
  };

  template <bool anything>
  struct lagrangian_helper<true, anything>
  {
    static std::vector<DomainType> get(const Implementation& imp)
    {
      const auto lps = imp.lagrangePoints();
      std::vector<DomainType> ret(lps.size());
      for (size_t ii = 0; ii < lps.size(); ++ii)
        ret[ii] = lps[ii].point();
      return ret;
    }
  };

public:
  static std::vector<DomainType> get(const Implementation& imp)
  {
    return lagrangian_helper<>::get(imp);
  }
}; // class LocalFiniteElementInterpolationLagrangepointsAccessor


} // namespace internal


template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementBasisWrapper : public LocalFiniteElementBasisInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementBasisWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  LocalFiniteElementBasisWrapper(Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
  {
  }

  LocalFiniteElementBasisWrapper(const Implementation& imp)
    : imp_(imp)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementBasisWrapper(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
  {
  }

  int order() const override final
  {
    return XT::Common::numeric_cast<int>(imp_.access().order());
  }

  size_t size() const override final
  {
    return XT::Common::numeric_cast<size_t>(imp_.access().size());
  }

  std::vector<RangeType> evaluate(const DomainType& xx) const override final
  {
    std::vector<RangeType> ret(size(), RangeType(0));
    imp_.access().evaluateFunction(xx, ret);
    return ret;
  }

  std::vector<JacobianRangeType> jacobian(const DomainType& xx) const override final
  {
    std::vector<JacobianRangeType> ret(size(), JacobianRangeType(0));
    imp_.access().evaluateJacobian(xx, ret);
    return ret;
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
}; // class LocalFiniteElementBasisWrapper


template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementInterpolationWrapper : public LocalFiniteElementInterpolationInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementInterpolationWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterpolationInterface<D, d, R, r, rC>;

  // This is what dune-localfunctions expects for interpolation.
  struct FunctionWrapper
  {
    // really?
    struct Traits
    {
      using DomainType = typename BaseType::DomainType;
      using RangeType = typename BaseType::RangeType;
    };
    // really!
    using DomainType = typename Traits::DomainType;
    using RangeType = typename Traits::RangeType;

    FunctionWrapper(const std::function<RangeType(DomainType)>& im)
      : imp(im)
    {
    }

    void evaluate(const DomainType& xx, RangeType& ret) const
    {
      ret = imp(xx);
    }

    const std::function<RangeType(DomainType)>& imp;
  };

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  LocalFiniteElementInterpolationWrapper(Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
  {
  }

  LocalFiniteElementInterpolationWrapper(const Implementation& imp)
    : imp_(imp)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementInterpolationWrapper(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
  {
  }

  std::vector<R> interpolate(const std::function<RangeType(DomainType)>& local_function) const override final
  {
    std::vector<R> ret;
    imp_.access().interpolate(FunctionWrapper(local_function), ret);
    return ret;
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
}; // class LocalFiniteElementInterpolationWrapper


template <class Implementation>
class LocalFiniteElementCoefficientsWrapper : public LocalFiniteElementCoefficientsInterface
{
  using ThisType = LocalFiniteElementCoefficientsWrapper;

public:
  LocalFiniteElementCoefficientsWrapper(Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
  {
  }

  LocalFiniteElementCoefficientsWrapper(const Implementation& imp)
    : imp_(imp)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementCoefficientsWrapper(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
  {
  }

  size_t size() const override final
  {
    return XT::Common::numeric_cast<size_t>(imp_.access().size());
  }

  const LocalKey& local_key(const size_t ii) const override final
  {
    if (ii >= size())
      DUNE_THROW(Exceptions::finite_element_error, "ii = " << ii << "\n   size() = " << size());
    return imp_.access().localKey(ii);
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
}; // class LocalFiniteElementCoefficientsWrapper


template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementWrapper : public LocalFiniteElementInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterface<D, d, R, r, rC>;

  using BasisWrapperType =
      LocalFiniteElementBasisWrapper<std::decay_t<typename Implementation::Traits::LocalBasisType>, D, d, R, r, rC>;
  using CoefficientsWrapperType =
      LocalFiniteElementCoefficientsWrapper<std::decay_t<typename Implementation::Traits::LocalCoefficientsType>>;
  using InterpolationWrapperType =
      LocalFiniteElementInterpolationWrapper<std::decay_t<typename Implementation::Traits::LocalInterpolationType>,
                                             D,
                                             d,
                                             R,
                                             r,
                                             rC>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::BasisType;
  using typename BaseType::CoefficientsType;
  using typename BaseType::InterpolationType;

private:
  using LpAccessor = internal::LocalFiniteElementInterpolationLagrangepointsAccessor<
      typename Implementation::Traits::LocalInterpolationType,
      DomainType>;

public:
  LocalFiniteElementWrapper(Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
  {
  }

  LocalFiniteElementWrapper(const Implementation& imp)
    : imp_(imp)
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
  {
  }

  template <class... Args>
  explicit LocalFiniteElementWrapper(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , basis_(imp_.access().localBasis())
    , coefficients_(imp_.access().localCoefficients())
    , interpolation_(imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
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
    return lagrange_points_.size() > 0;
  }

  const std::vector<DomainType>& lagrange_points() const override final
  {
    if (!is_lagrangian())
      DUNE_THROW(Exceptions::finite_element_error, "do not call lagrange_points() if is_lagrangian() is false!");
    return lagrange_points_;
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
  const BasisWrapperType basis_;
  const CoefficientsWrapperType coefficients_;
  const InterpolationWrapperType interpolation_;
  const std::vector<DomainType> lagrange_points_;
}; // class LocalFiniteElementWrapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH
