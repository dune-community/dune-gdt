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
  using typename BaseType::DerivativeRangeType;

  LocalFiniteElementBasisWrapper(const GeometryType& gt, Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
    , geometry_type_(gt)
  {
  }

  LocalFiniteElementBasisWrapper(const GeometryType& gt, const Implementation& imp)
    : imp_(imp)
    , geometry_type_(gt)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementBasisWrapper(const GeometryType& gt, Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , geometry_type_(gt)
  {
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  int order() const override final
  {
    return XT::Common::numeric_cast<int>(imp_.access().order());
  }

  size_t size() const override final
  {
    return XT::Common::numeric_cast<size_t>(imp_.access().size());
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& point_in_reference_element, std::vector<RangeType>& result) const override final
  {
    imp_.access().evaluateFunction(point_in_reference_element, result);
  }

  using BaseType::jacobian;

  void jacobian(const DomainType& point_in_reference_element,
                std::vector<DerivativeRangeType>& result) const override final
  {
    imp_.access().evaluateJacobian(point_in_reference_element, result);
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
  const GeometryType geometry_type_;
}; // class LocalFiniteElementBasisWrapper


template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementInterpolationWrapper : public LocalFiniteElementInterpolationInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementInterpolationWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterpolationInterface<D, d, R, r, rC>;

  template <class R_ = R, size_t r_ = r, size_t rC_ = rC>
  struct RangeTypeSelector
  {
    using type = FieldMatrix<R_, r_, rC_>;
  };

  template <class R_, size_t r_>
  struct RangeTypeSelector<R_, r_, 1>
  {
    using type = FieldVector<R_, r_>;
  };

  // This is what dune-localfunctions expects for interpolation.
  struct FunctionWrapper
  {
    // really?
    struct Traits
    {
      using DomainType = typename BaseType::DomainType;
      using RangeType = typename RangeTypeSelector<>::type;
    };
    // really!
    using DomainType = typename Traits::DomainType;
    using RangeType = typename RangeTypeSelector<>::type;

    FunctionWrapper(const std::function<RangeType(DomainType)>& im)
      : imp(im)
    {
    }

    void evaluate(const DomainType& point_in_reference_element, RangeType& ret) const
    {
      ret = imp(point_in_reference_element);
    }

    const std::function<RangeType(DomainType)>& imp;
  };

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  LocalFiniteElementInterpolationWrapper(const GeometryType& gt, Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
    , geometry_type_(gt)
  {
  }

  LocalFiniteElementInterpolationWrapper(const GeometryType& gt, const Implementation& imp)
    : imp_(imp)
    , geometry_type_(gt)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementInterpolationWrapper(const GeometryType& gt, Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , geometry_type_(gt)
  {
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  using BaseType::interpolate;

  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   std::vector<R>& result) const override final
  {
    imp_.access().interpolate(FunctionWrapper(local_function), result);
  }

private:
  const XT::Common::ConstStorageProvider<Implementation> imp_;
  const GeometryType geometry_type_;
}; // class LocalFiniteElementInterpolationWrapper


template <class Implementation, class D, size_t d>
class LocalFiniteElementCoefficientsWrapper : public LocalFiniteElementCoefficientsInterface<D, d>
{
  using ThisType = LocalFiniteElementCoefficientsWrapper<Implementation, D, d>;

public:
  LocalFiniteElementCoefficientsWrapper(const GeometryType& gt, Implementation*&& imp_ptr)
    : imp_(std::move(imp_ptr))
    , geometry_type_(gt)
  {
  }

  LocalFiniteElementCoefficientsWrapper(const GeometryType& gt, const Implementation& imp)
    : imp_(imp)
    , geometry_type_(gt)
  {
  }

  template <class... Args>
  explicit LocalFiniteElementCoefficientsWrapper(const GeometryType& gt, Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , geometry_type_(gt)
  {
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
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
  const GeometryType geometry_type_;
}; // class LocalFiniteElementCoefficientsWrapper


template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementWrapper : public LocalFiniteElementInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterface<D, d, R, r, rC>;

  using BasisWrapperType =
      LocalFiniteElementBasisWrapper<std::decay_t<typename Implementation::Traits::LocalBasisType>, D, d, R, r, rC>;
  using CoefficientsWrapperType =
      LocalFiniteElementCoefficientsWrapper<std::decay_t<typename Implementation::Traits::LocalCoefficientsType>, D, d>;
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
    , geometry_type_(imp_.access().type())
    , basis_(imp_.access().type(), imp_.access().localBasis())
    , coefficients_(imp_.access().type(), imp_.access().localCoefficients())
    , interpolation_(imp_.access().type(), imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
  {
  }

  LocalFiniteElementWrapper(const Implementation& imp)
    : imp_(imp)
    , geometry_type_(imp_.access().type())
    , basis_(imp_.access().type(), imp_.access().localBasis())
    , coefficients_(imp_.access().type(), imp_.access().localCoefficients())
    , interpolation_(imp_.access().type(), imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
  {
  }

  template <class... Args>
  explicit LocalFiniteElementWrapper(Args&&... args)
    : imp_(new Implementation(std::forward<Args>(args)...))
    , geometry_type_(imp_.access().type())
    , basis_(imp_.access().type(), imp_.access().localBasis())
    , coefficients_(imp_.access().type(), imp_.access().localCoefficients())
    , interpolation_(imp_.access().type(), imp_.access().localInterpolation())
    , lagrange_points_(LpAccessor::get(imp_.access().localInterpolation()))
  {
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
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
  const GeometryType geometry_type_;
  const BasisWrapperType basis_;
  const CoefficientsWrapperType coefficients_;
  const InterpolationWrapperType interpolation_;
  const std::vector<DomainType> lagrange_points_;
}; // class LocalFiniteElementWrapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH
