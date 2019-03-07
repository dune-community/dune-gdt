// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH

#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/type_traits.hh>

#include "default.hh"
#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace internal {


/**
 * Can be used by wrappers of local finite elements from dune-localfunctions to obtain lagrange points, if those exist.
 *
 * Most likely, you do not want to use this class directly.
 */
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
      for (unsigned int ii = 0; ii < lps.size(); ++ii)
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


/**
 * \brief Wraps a local finite element from dune-localfunctions to implement LocalFiniteElementBasisInterface.
 *
 * Most likely, you do not want to use this class directly, but LocalFiniteElementWrapper.
 *
 * \sa   LocalFiniteElementBasisInterface
 * \sa   LocalFiniteElementWrapper
 */
template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementBasisWrapper : public LocalFiniteElementBasisInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementBasisWrapper<Implementation, D, d, R, r, rC>;
  using BaseType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;

public:
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  LocalFiniteElementBasisWrapper(std::shared_ptr<const Implementation> imp)
    : imp_(imp)
    , geometry_type_(imp_->type())
  {}

  LocalFiniteElementBasisWrapper(Implementation*&& imp_ptr)
    : LocalFiniteElementBasisWrapper(std::shared_ptr<Implementation>(std::move(imp_ptr)))
  {}

  template <class ImpType,
            typename = typename std::enable_if<std::is_same<ImpType, Implementation>::value
                                               && std::is_copy_constructible<ImpType>::value>::type>
  LocalFiniteElementBasisWrapper(const ImpType& imp)
    : LocalFiniteElementBasisWrapper(new ImpType(imp))
  {}

  template <class... Args>
  explicit LocalFiniteElementBasisWrapper(Args&&... args)
    : LocalFiniteElementBasisWrapper(new Implementation(std::forward<Args>(args)...))
  {}

  LocalFiniteElementBasisWrapper(const ThisType& other) = default;

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
    return XT::Common::numeric_cast<int>(imp_->localBasis().order());
  }

  size_t size() const override final
  {
    return XT::Common::numeric_cast<size_t>(imp_->localBasis().size());
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& point_in_reference_element, std::vector<RangeType>& result) const override final
  {
    imp_->localBasis().evaluateFunction(point_in_reference_element, result);
  }

  using BaseType::jacobian;

  void jacobian(const DomainType& point_in_reference_element,
                std::vector<DerivativeRangeType>& result) const override final
  {
    imp_->localBasis().evaluateJacobian(point_in_reference_element, result);
  }

private:
  const std::shared_ptr<const Implementation> imp_;
  const GeometryType geometry_type_;
}; // class LocalFiniteElementBasisWrapper


/**
 * \brief Wraps a local finite element from dune-localfunctions to implement LocalFiniteElementInterpolationInterface.
 *
 * Most likely, you do not want to use this class directly, but LocalFiniteElementWrapper.
 *
 * \sa   LocalFiniteElementInterpolationInterface
 * \sa   LocalFiniteElementWrapper
 */
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
    {}

    void evaluate(const DomainType& point_in_reference_element, RangeType& ret) const
    {
      ret = imp(point_in_reference_element);
    }

    const std::function<RangeType(DomainType)>& imp;
  };

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  LocalFiniteElementInterpolationWrapper(std::shared_ptr<const Implementation> imp)
    : imp_(imp)
    , geometry_type_(imp_->type())
    , dofs_(imp_->size())
  {}

  LocalFiniteElementInterpolationWrapper(Implementation*&& imp_ptr)
    : LocalFiniteElementInterpolationWrapper(std::shared_ptr<Implementation>(std::move(imp_ptr)))
  {}

  template <class ImpType,
            typename = typename std::enable_if<std::is_same<ImpType, Implementation>::value
                                               && std::is_copy_constructible<ImpType>::value>::type>
  LocalFiniteElementInterpolationWrapper(const ImpType& imp)
    : LocalFiniteElementInterpolationWrapper(new ImpType(imp))
  {}

  template <class... Args>
  explicit LocalFiniteElementInterpolationWrapper(Args&&... args)
    : LocalFiniteElementInterpolationWrapper(new Implementation(std::forward<Args>(args)...))
  {}

  LocalFiniteElementInterpolationWrapper(const ThisType& other) = default;

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
    return XT::Common::numeric_cast<size_t>(imp_->size());
  }

  using BaseType::interpolate;

  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   const int /*order*/,
                   DynamicVector<R>& dofs) const override final
  {
    imp_->localInterpolation().interpolate(FunctionWrapper(local_function), dofs_);
    const size_t sz = this->size();
    if (dofs.size() != sz)
      dofs.resize(sz);
    for (size_t ii = 0; ii < sz; ++ii)
      dofs[ii] = dofs_[ii];
  }

private:
  const std::shared_ptr<const Implementation> imp_;
  const GeometryType geometry_type_;
  mutable std::vector<R> dofs_;
}; // class LocalFiniteElementInterpolationWrapper


/**
 * \brief Wraps a local finite element from dune-localfunctions to implement LocalFiniteElementCoefficientsInterface.
 *
 * Most likely, you do not want to use this class directly, but LocalFiniteElementWrapper.
 *
 * \sa   LocalFiniteElementCoefficientsInterface
 * \sa   LocalFiniteElementWrapper
 */
template <class Implementation, class D, size_t d>
class LocalFiniteElementCoefficientsWrapper : public LocalFiniteElementCoefficientsInterface<D, d>
{
  using ThisType = LocalFiniteElementCoefficientsWrapper<Implementation, D, d>;
  using BaseType = LocalFiniteElementCoefficientsInterface<D, d>;

public:
  LocalFiniteElementCoefficientsWrapper(std::shared_ptr<const Implementation> imp)
    : imp_(imp)
    , geometry_type_(imp_->type())
  {}

  LocalFiniteElementCoefficientsWrapper(Implementation*&& imp_ptr)
    : LocalFiniteElementCoefficientsWrapper(std::shared_ptr<Implementation>(std::move(imp_ptr)))
  {}

  template <class ImpType,
            typename = typename std::enable_if<std::is_same<ImpType, Implementation>::value
                                               && std::is_copy_constructible<ImpType>::value>::type>
  LocalFiniteElementCoefficientsWrapper(const ImpType& imp)
    : LocalFiniteElementCoefficientsWrapper(new ImpType(imp))
  {}

  template <class... Args>
  explicit LocalFiniteElementCoefficientsWrapper(Args&&... args)
    : LocalFiniteElementCoefficientsWrapper(new Implementation(std::forward<Args>(args)...))
  {}

  LocalFiniteElementCoefficientsWrapper(const ThisType& other) = default;

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
    return XT::Common::numeric_cast<size_t>(imp_->localCoefficients().size());
  }

  const LocalKey& local_key(const size_t ii) const override final
  {
    if (ii >= size())
      DUNE_THROW(Exceptions::finite_element_error, "ii = " << ii << "\n   size() = " << size());
    return imp_->localCoefficients().localKey(XT::Common::numeric_cast<int>(ii));
  }

private:
  const std::shared_ptr<const Implementation> imp_;
  const GeometryType geometry_type_;
}; // class LocalFiniteElementCoefficientsWrapper


/**
 * \brief     Wraps a local finite element from dune-localfunctions.
 *
 * \attention We presume that the Implementation from dune-localfunctions is thread-safe (for virtual copyability)!
 *
 * If thread-safety cannot be ensured, we will need to think about generators (since some of the ingredients are not
 * copyable).
 *
 * \sa        LocalFiniteElementBasisWrapper
 * \sa        LocalFiniteElementInterpolationWrapper
 * \sa        LocalFiniteElementCoefficientsWrapper
 * \sa        LocalFiniteElementInterface
 */
template <class Implementation, class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementWrapper
  : XT::Common::ConstSharedStorageProvider<Implementation>
  , public LocalFiniteElementDefault<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementWrapper<Implementation, D, d, R, r, rC>;
  using ImplementationProvider = XT::Common::ConstSharedStorageProvider<Implementation>;
  using BaseType = LocalFiniteElementDefault<D, d, R, r, rC>;

public:
  using typename BaseType::BasisType;
  using typename BaseType::CoefficientsType;
  using typename BaseType::DomainType;
  using typename BaseType::InterpolationType;

private:
  using LpAccessor = internal::LocalFiniteElementInterpolationLagrangepointsAccessor<
      typename Implementation::Traits::LocalInterpolationType,
      DomainType>;

public:
  LocalFiniteElementWrapper(const int ord, Implementation*&& imp_ptr)
    : ImplementationProvider(std::move(imp_ptr))
    , BaseType(
          ord,
          new LocalFiniteElementBasisWrapper<Implementation, D, d, R, r, rC>(ImplementationProvider::access()),
          new LocalFiniteElementCoefficientsWrapper<Implementation, D, d>(ImplementationProvider::access()),
          new LocalFiniteElementInterpolationWrapper<Implementation, D, d, R, r, rC>(ImplementationProvider::access()),
          LpAccessor::get(ImplementationProvider::access()->localInterpolation()))
  {}

  template <class ImpType,
            typename = typename std::enable_if<std::is_same<ImpType, Implementation>::value
                                               && std::is_copy_constructible<ImpType>::value>::type>
  LocalFiniteElementWrapper(const int ord, const ImpType& imp)
    : LocalFiniteElementWrapper(ord, new ImpType(imp))
  {}

  template <class... Args>
  explicit LocalFiniteElementWrapper(const int ord, Args&&... args)
    : LocalFiniteElementWrapper(ord, new Implementation(std::forward<Args>(args)...))
  {}
}; // class LocalFiniteElementWrapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_WRAPPER_HH
