// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_FLATTOP_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_FLATTOP_HH

#include <functional>

#include <dune/gdt/exceptions.hh>

#include "default.hh"
#include "interfaces.hh"
#include "lagrange.hh"
#include "power.hh"

namespace Dune {
namespace GDT {


/**
 * Inspired by [Brenner, Davis, Sung, 2014, A partition of unity method for the displacement obstacle problem of clamped
 *              Kirchhoff plates], section 2
 *
 * \sa LocalFlatTop2dCubeFiniteElement
 * \sa LocalFlatTopFiniteElementFactory
 */
template <class D = double, class R = double>
class LocalFlatTop2dCubeFiniteElementBasis : public LocalFiniteElementBasisInterface<D, 2, R, 1, 1>
{
  using ThisType = LocalFlatTop2dCubeFiniteElementBasis<D, R>;
  using BaseType = LocalFiniteElementBasisInterface<D, 2, R, 1, 1>;

public:
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  LocalFlatTop2dCubeFiniteElementBasis(const double& overlap = 0.5)
    : geometry_type_(GeometryType::cube, 2)
  {
    DUNE_THROW_IF(
        !(overlap > 0.), Exceptions::finite_element_error, "Overlap has to be in (0, 1], is " << overlap << "!");
    DUNE_THROW_IF(overlap > 1., Exceptions::finite_element_error, "Overlap has to be in (0, 1], is " << overlap << "!");
    // we cannot let L_ and R_ be members and define phi_L_ and phi_R_ in the ctor initializer list, as they will
    // copy/reference broken/empty/default L_ and R_
    const auto L_ = (1. - overlap) / 2.;
    const auto R_ = (1. + overlap) / 2.;
    phi_L_ = [=](const D& x) {
      if (x < L_)
        return 1.;
      else if (x > R_)
        return 0.;
      else
        return -1. * (x / overlap) + R_ / overlap;
    };
    grad_phi_L_ = [=](const D& x) {
      if (x < L_)
        return 0.;
      else if (x > R_)
        return 0.;
      else
        return -1. / overlap;
    };
    phi_R_ = [=](const D& x) {
      if (x < L_)
        return 0.;
      else if (x > R_)
        return 1.;
      else
        return (x / overlap) - L_ / overlap;
    };
    grad_phi_R_ = [=](const D& x) {
      if (x < L_)
        return 0.;
      else if (x > R_)
        return 0.;
      else
        return (1. / overlap);
    };
  } // LocalFlatTop2dCubeFiniteElementBasis(...)

  LocalFlatTop2dCubeFiniteElementBasis(const ThisType& other) = default;
  LocalFlatTop2dCubeFiniteElementBasis(ThisType&& source) = default;

  ThisType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return geometry_type_;
  }

  int order() const override final
  {
    return 2;
  }

  size_t size() const override final
  {
    return 4;
  }

  using BaseType::evaluate;

  void evaluate(const DomainType& point_in_reference_element, std::vector<RangeType>& result) const override final
  {
    if (result.size() < 4)
      result.resize(4);
    const auto& x = point_in_reference_element[0];
    const auto& y = point_in_reference_element[1];
    // each shape function is associated with one corner, the positions of which are known on the reference element
    // * shape function 0: bottom left corner
    result[0] = phi_L_(x) * phi_L_(y);
    // * shape function 1: bottom right corner
    result[1] = phi_R_(x) * phi_L_(y);
    // * shape function 2: top left corner
    result[2] = phi_L_(x) * phi_R_(y);
    // * shape function 3: top right corner
    result[3] = phi_R_(x) * phi_R_(y);
  } // ... evaluate(...)

  using BaseType::jacobian;

  void jacobian(const DomainType& point_in_reference_element,
                std::vector<DerivativeRangeType>& result) const override final
  {
    if (result.size() < 4)
      result.resize(4);
    const auto& x = point_in_reference_element[0];
    const auto& y = point_in_reference_element[1];
    // * shape function 0
    result[0][0][0] = grad_phi_L_(x) * phi_L_(y);
    result[0][0][1] = phi_L_(x) * grad_phi_L_(y);
    // * shape function 1
    result[1][0][0] = grad_phi_R_(x) * phi_L_(y);
    result[1][0][1] = phi_R_(x) * grad_phi_L_(y);
    // * shape function 2
    result[2][0][0] = grad_phi_L_(x) * phi_R_(y);
    result[2][0][1] = phi_L_(x) * grad_phi_R_(y);
    // * shape function 3
    result[3][0][0] = grad_phi_R_(x) * phi_R_(y);
    result[3][0][1] = phi_R_(x) * grad_phi_R_(y);
  } // ... jacobian(...)

private:
  const GeometryType geometry_type_;
  std::function<R(const D&)> phi_L_;
  std::function<R(const D&)> phi_R_;
  std::function<R(const D&)> grad_phi_L_;
  std::function<R(const D&)> grad_phi_R_;
}; // class LocalFlatTop2dCubeFiniteElementBasis


template <class D = double, class R = double>
class LocalFlatTop2dCubeFiniteElement : public LocalFiniteElementDefault<D, 2, R, 1>
{
  using ThisType = LocalFlatTop2dCubeFiniteElement<D, R>;
  using BaseType = LocalFiniteElementDefault<D, 2, R, 1>;

public:
  LocalFlatTop2dCubeFiniteElement(const double& overlap = 0.5)
    : BaseType(
          1,
          LocalFlatTop2dCubeFiniteElementBasis<D, R>(overlap).copy(),
          LocalLagrangeFiniteElementFactory<D, 2, R, 1>::create(GeometryType(GeometryType::cube, 2), 1)
              ->coefficients()
              .copy(),
          LocalL2FiniteElementInterpolation<D, 2, R, 1>(LocalFlatTop2dCubeFiniteElementBasis<D, R>(overlap)).copy(),
          {})
  {}
}; // class LocalFlatTop2dCubeFiniteElement


template <class D, size_t d, class R, size_t r = 1>
class LocalFlatTopFiniteElementFactory
{
  using ScalarLocalFiniteElementType = LocalFiniteElementInterface<D, d, R, 1>;

public:
  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, r>;

private:
  static std::string order_error(const GeometryType& geometry_type, const int order)
  {
    std::stringstream ss;
    ss << "when creating a local FlatTop finite element: the FlatTopLocalFiniteElement is known to fail in " << d
       << "d on a " << geometry_type << " reference element for order " << order
       << " (if you think it is working, update this check)!";
    return ss.str();
  }

  static std::string geometry_error(const GeometryType& geometry_type, const int order)
  {
    std::stringstream ss;
    ss << "when creating a local FlatTop finite element: this is untested!\n"
       << "Please update this check if you believe that FlatTopLocalFiniteElement is available for\n- dimension: " << d
       << "\n- geometry_type: " << geometry_type << "\n- order: " << order;
    return ss.str();
  }

  // Fist we create the scalar FE ...

  template <size_t d_ = d, bool anything = true>
  struct scalar_helper
  {
    static std::unique_ptr<ScalarLocalFiniteElementType>
    create(const GeometryType& geometry_type, const int& order, const D& overlap = 0.5)
    {
      DUNE_THROW_IF(geometry_type.dim() != d,
                    Exceptions::finite_element_error,
                    "geometry_type = " << geometry_type << "\nd = " << d);
      // checks
      if (d == 2) {
        if (geometry_type == GeometryTypes::cube(2))
          DUNE_THROW_IF(order != 1, Exceptions::finite_element_error, order_error(geometry_type, order));
        else
          DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      } else
        DUNE_THROW(Exceptions::finite_element_error, geometry_error(geometry_type, order));
      // the actual finite element
      return std::unique_ptr<LocalFiniteElementInterface<D, d, R, 1>>(
          new LocalFlatTop2dCubeFiniteElement<D, R>(overlap));
    }
  }; // helper<...>

  template <bool anything>
  struct scalar_helper<0, anything>
  {
    static std::unique_ptr<ScalarLocalFiniteElementType>
    create(const GeometryType& geometry_type, const int& /*order*/, const D& /*overlap*/ = 0.5)
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
    static std::unique_ptr<LocalFiniteElementType>
    create(const GeometryType& geometry_type, const int& order, const D& overlap = 0.5)
    {
      return make_local_powered_finite_element<r>(*scalar_helper<>::create(geometry_type, order, overlap));
    }
  };

  template <size_t d_>
  struct helper<d_, 1>
  {
    static std::unique_ptr<LocalFiniteElementType>
    create(const GeometryType& geometry_type, const int& order, const D& overlap = 0.5)
    {
      return scalar_helper<>::create(geometry_type, order, overlap);
    }
  };

public:
  static std::unique_ptr<LocalFiniteElementInterface<D, d, R, r>>
  create(const GeometryType& geometry_type, const int& order, const D& overlap = 0.5)
  {
    return helper<>::create(geometry_type, order, overlap);
  }
}; // class LocalFlatTopFiniteElementFactory


template <class D, size_t d, class R, size_t r = 1>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, r>>
make_local_flattop_finite_element(const GeometryType& geometry_type, const int order, const D& overlap = 0.5)
{
  return LocalFlatTopFiniteElementFactory<D, d, R, r>::create(geometry_type, order, overlap);
}


template <class D, size_t d, class R, size_t r = 1>
class LocalFlatTopFiniteElementFamily : public ThreadSafeDefaultLocalFiniteElementFamily<D, d, R, r>
{
  using BaseType = ThreadSafeDefaultLocalFiniteElementFamily<D, d, R, r>;

public:
  LocalFlatTopFiniteElementFamily(const D& overlap = 0.5)
    : BaseType([=](const auto& geometry_type, const auto& order) {
      return LocalFlatTopFiniteElementFactory<D, d, R, r>::create(geometry_type, order, overlap);
    })
  {}
}; // ... LocalFlatTopFiniteElementFamily(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_FLATTOP_HH
