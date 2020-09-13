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

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH

#include <sstream>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/raviartthomas.hh>

#include <dune/xt/la/container/common/matrix/dense.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/la/solver.hh>

#include <dune/gdt/exceptions.hh>

#include "default.hh"
#include "interfaces.hh"
#include "orthonormal.hh"
#include "wrapper.hh"
#include "0d.hh"

namespace Dune {
namespace GDT {


/**
 * The one from dune-localfunction is broken for some reference elements, at least for 2d and 3d cubes.
 */
template <class D, size_t d, class R>
class LocalRaviartThomasInterpolation : public LocalFiniteElementInterpolationInterface<D, d, R, d, 1>
{
  using ThisType = LocalRaviartThomasInterpolation;
  using BaseType = LocalFiniteElementInterpolationInterface<D, d, R, d, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  using BasisType = LocalFiniteElementBasisInterface<D, d, R, d>;
  using CoefficientsType = LocalFiniteElementCoefficientsInterface<D, d>;

  LocalRaviartThomasInterpolation(const int ord, const BasisType& basis, const CoefficientsType& coefficients)
    : order_(ord)
    , basis_(basis.copy())
    , coefficients_(coefficients.copy())
  {
    if (basis_->size() != coefficients_->size())
      DUNE_THROW(Exceptions::finite_element_error,
                 "basis_->size() = " << basis_->size() << "\n   "
                                     << "coefficients_->size() = " << coefficients_->size());
  }

  LocalRaviartThomasInterpolation(const ThisType& other)
    : order_(other.order_)
    , basis_(other.basis_->copy())
    , coefficients_(other.coefficients_->copy())
  {}

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  const GeometryType& geometry_type() const override final
  {
    return basis_->geometry_type();
  }

  size_t size() const override final
  {
    return basis_->size();
  }

  using BaseType::interpolate;

  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   const int order,
                   DynamicVector<R>& dofs) const override final
  {
    // prepare
    const size_t sz = this->size();
    if (dofs.size() < sz)
      dofs.resize(sz);
    std::vector<bool> local_key_was_handled(sz, false);
    const auto& reference_element = ReferenceElements<D, d>::general(this->geometry_type());
    // determine the face dofs, therefore walk the intersections
    const auto intersection_to_local_key_map = coefficients_->local_key_indices(1);
    for (int intersection_index = 0; intersection_index < reference_element.size(1); ++intersection_index) {
      const auto& intersection_geometry_type = reference_element.type(intersection_index, 1);
      const auto& local_keys_assosiated_with_intersection = intersection_to_local_key_map[intersection_index];
      if (local_keys_assosiated_with_intersection.size() > 0) {
        const auto intersection_fe =
            make_local_orthonormal_finite_element<D, d - 1, R>(intersection_geometry_type, order_);
        const auto& intersection_Pk_basis = intersection_fe->basis();
        DUNE_THROW_IF(intersection_Pk_basis.size() != local_keys_assosiated_with_intersection.size(),
                      Exceptions::finite_element_error,
                      error_msg_prefix() << "intersection_Pk_basis.size() = " << intersection_Pk_basis.size()
                                         << "\n   local_keys_assosiated_with_intersection.size() = "
                                         << local_keys_assosiated_with_intersection.size());
        XT::LA::CommonDenseMatrix<R> lhs(
            local_keys_assosiated_with_intersection.size(), intersection_Pk_basis.size(), 0);
        XT::LA::CommonDenseVector<R> rhs(intersection_Pk_basis.size(), 0);
        // do a face quadrature
        for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                 intersection_geometry_type,
                 std::max(basis_->order() + intersection_Pk_basis.order(), order + intersection_Pk_basis.order()))) {
          const auto point_on_reference_intersection = quadrature_point.position();
          const auto point_in_reference_element =
              reference_element.template geometry<1>(intersection_index).global(point_on_reference_intersection);
          const auto quadrature_weight = quadrature_point.weight();
          const auto integration_outer_normal = reference_element.integrationOuterNormal(intersection_index);
          const auto rt_basis_values = basis_->evaluate(point_in_reference_element);
          const auto intersection_Pk_basis_values = intersection_Pk_basis.evaluate(point_on_reference_intersection);
          const auto local_function_values = local_function(point_in_reference_element);
          for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
            const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
            for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
              lhs.add_to_entry(ii,
                               jj,
                               quadrature_weight * (rt_basis_values[local_key_index] * integration_outer_normal)
                                   * intersection_Pk_basis_values[jj]);
          }
          for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
            rhs[jj] += quadrature_weight * (local_function_values * integration_outer_normal)
                       * intersection_Pk_basis_values[jj];
        }
        XT::LA::CommonDenseVector<R> intersection_dofs(local_keys_assosiated_with_intersection.size(), 0);
        try {
          intersection_dofs = XT::LA::solve(lhs, rhs);
        } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
          DUNE_THROW(Exceptions::finite_element_error,
                     error_msg_prefix() << "Failed to solve for DoFs associated with intersection "
                                        << intersection_index << ", this was the original error:\n   " << ee.what());
        }
        for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
          const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
          assert(!local_key_was_handled[local_key_index]);
          dofs[local_key_index] = intersection_dofs[ii];
          local_key_was_handled[local_key_index] = true;
        }
      }
    }
    // determine the volume dofs
    const auto element_to_local_key_map = coefficients_->local_key_indices(0);
    const auto& local_keys_assosiated_with_element = element_to_local_key_map[0];
    if (local_keys_assosiated_with_element.size() > 0) {
      DUNE_THROW_IF(order_ < 1,
                    Exceptions::finite_element_error,
                    error_msg_prefix() << "DoFs associated with the element only make sense for orders >= 1!");
      const auto element_fe = make_local_orthonormal_finite_element<D, d, R, d>(this->geometry_type(), order_ - 1);
      const auto& element_Pkminus1_basis = element_fe->basis();
      DUNE_THROW_IF(element_Pkminus1_basis.size() != local_keys_assosiated_with_element.size(),
                    Exceptions::finite_element_error,
                    error_msg_prefix() << "element_Pkminus1_basis.size() = " << element_Pkminus1_basis.size()
                                       << "\n   local_keys_assosiated_with_element.size() = "
                                       << local_keys_assosiated_with_element.size());
      XT::LA::CommonDenseMatrix<R> lhs(local_keys_assosiated_with_element.size(), element_Pkminus1_basis.size(), 0);
      XT::LA::CommonDenseVector<R> rhs(element_Pkminus1_basis.size(), 0);
      // do a volume quadrature
      for (auto&& quadrature_point : QuadratureRules<D, d>::rule(
               this->geometry_type(),
               std::max(basis_->order() + element_Pkminus1_basis.order(), order + element_Pkminus1_basis.order()))) {
        const auto point_in_reference_element = quadrature_point.position();
        const auto quadrature_weight = quadrature_point.weight();
        const auto rt_basis_values = basis_->evaluate(point_in_reference_element);
        const auto element_Pkminus1_basis_values = element_Pkminus1_basis.evaluate(point_in_reference_element);
        const auto local_function_values = local_function(point_in_reference_element);
        for (size_t ii = 0; ii < local_keys_assosiated_with_element.size(); ++ii) {
          const size_t local_key_index = local_keys_assosiated_with_element[ii];
          for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
            lhs.add_to_entry(
                ii, jj, quadrature_weight * (rt_basis_values[local_key_index] * element_Pkminus1_basis_values[jj]));
        }
        for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
          rhs[jj] += quadrature_weight * (local_function_values * element_Pkminus1_basis_values[jj]);
      }
      XT::LA::CommonDenseVector<R> element_dofs(local_keys_assosiated_with_element.size(), 0);
      try {
        element_dofs = XT::LA::solve(lhs, rhs);
      } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
        DUNE_THROW(Exceptions::finite_element_error,
                   error_msg_prefix() << "Failed to solve for volume DoFs, this was the original error:\n   "
                                      << ee.what());
      }
      for (size_t ii = 0; ii < local_keys_assosiated_with_element.size(); ++ii) {
        const size_t local_key_index = local_keys_assosiated_with_element[ii];
        assert(!local_key_was_handled[local_key_index]);
        dofs[local_key_index] = element_dofs[ii];
        local_key_was_handled[local_key_index] = true;
      }
    }
    // final checks that there are no other dofs left
    for (size_t ii = 0; ii < sz; ++ii)
      DUNE_THROW_IF(
          !local_key_was_handled[ii],
          Exceptions::finite_element_error,
          error_msg_prefix() << "The following DoF is neither associated with an intersection, nor with the element!"
                             << "\n   local DoF index: " << ii
                             << "\n   associated local_key: " << coefficients_->local_key(ii));
  } // ... interpolate(...)

private:
  std::string error_msg_prefix() const
  {
    std::stringstream ss;
    ss << "when interpolating with a local Raviart-Thomas Finite Element of order " << order_ << " on a "
       << this->geometry_type() << " reference element!\n   ";
    return ss.str();
  }

  const int order_;
  const std::unique_ptr<BasisType> basis_;
  const std::unique_ptr<CoefficientsType> coefficients_;
}; // class LocalRaviartThomasInterpolation


/**
 * \note Update this class if anything changes in dune-localfunctions.
 */
template <class D, size_t d, class R>
class LocalRaviartThomasFiniteElementFactory
{
  static_assert(0 <= d && d <= 3, "There is no local Raviart-Thomas finite element available for other dimension!");

public:
  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, d, 1>;

private:
  template <size_t d_ = d, bool anything = true>
  struct helper;

  template <bool anything>
  struct helper<0, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& /*order*/)
    {
      // If we need this, and geometry_type.dim() == 0, we must simply implement the corresponding ctors of the 0d FE!
      DUNE_THROW_IF(
          geometry_type.dim() != 0 || !geometry_type.isSimplex(),
          Exceptions::finite_element_error,
          "when creating a local 0d orthonomal finite element: not available for geometry_type = " << geometry_type);
      return std::make_unique<Local0dFiniteElement<D, R>>();
    }
  };

  template <bool anything>
  struct helper<1, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      // everything is a simplex in 1d
      using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
      return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(order, new FE(geometry_type, order));
    }
  };

  template <bool anything>
  struct helper<2, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      if (geometry_type.isSimplex()) {
        using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
        return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(order, new FE(geometry_type, order));
      } else if (geometry_type.isCube()) {
        if (order == 0) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 0>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(0, new FE());
        } else if (order == 1) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 1>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(1, new FE());
        } else if (order == 2) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 2>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(2, new FE());
        } else if (order == 3) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 3>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(3, new FE());
        } else if (order == 4) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 4>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(4, new FE());
        } else
          DUNE_THROW(Exceptions::finite_element_error,
                     "when creating a local Raviart-Thomas finite element: there is none available on cubes for order "
                         << order << " (if you think there is, update this class)!");
      } else
        DUNE_THROW(Exceptions::finite_element_error,
                   "when creating a local Raviart-Thomas finite element: there is none available for the "
                   "following geometry type: "
                       << geometry_type << "(if you think there is, update this class)!");
    }
  }; // helper<2, ...>

  template <bool anything>
  struct helper<3, anything>
  {
    static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
    {
      if (geometry_type.isSimplex()) {
        using FE = RaviartThomasSimplexLocalFiniteElement<d, D, R>;
        return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(order, new FE(geometry_type, order));
      } else if (geometry_type.isCube()) {
        if (order == 0) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 0>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(0, new FE());
        } else if (order == 1) {
          using FE = RaviartThomasCubeLocalFiniteElement<D, R, d, 1>;
          return std::make_unique<LocalFiniteElementWrapper<FE, D, d, R, d>>(1, new FE());
        } else
          DUNE_THROW(Exceptions::finite_element_error,
                     "when creating a local Raviart-Thomas finite element: there is none available on cubes for order "
                         << order << " (if you think there is, update this class)!");
      } else
        DUNE_THROW(Exceptions::finite_element_error,
                   "when creating a local Raviart-Thomas finite element: there is none available for the "
                   "following geometry type: "
                       << geometry_type << "(if you think there is, update this class)!");
    }
  }; // helper<3, ...>

public:
  static std::unique_ptr<LocalFiniteElementType> create(const GeometryType& geometry_type, const int& order)
  {
    // create the finite element
    auto fe = helper<>::create(geometry_type, order);
    // but use our interpolation
    return std::make_unique<LocalFiniteElementDefault<D, d, R, d, 1>>(
        order,
        fe->basis().copy(),
        fe->coefficients().copy(),
        new LocalRaviartThomasInterpolation<D, d, R>(order, fe->basis(), fe->coefficients()));
  }
}; // class LocalRaviartThomasFiniteElementFactory


template <class D, size_t d, class R>
std::unique_ptr<LocalFiniteElementInterface<D, d, R, d, 1>>
make_local_raviart_thomas_finite_element(const GeometryType& geometry_type, const int order)
{
  return LocalRaviartThomasFiniteElementFactory<D, d, R>::create(geometry_type, order);
}


template <class D, size_t d, class R>
class LocalRaviartThomasFiniteElementFamily : public ThreadSafeDefaultLocalFiniteElementFamily<D, d, R, d, 1>
{
  using BaseType = ThreadSafeDefaultLocalFiniteElementFamily<D, d, R, d, 1>;

public:
  LocalRaviartThomasFiniteElementFamily()
    : BaseType([](const auto& geometry_type, const auto& order) {
      // Can't figure out why this lock is needed, but for some reasons without it we are not thread-safe
      static std::mutex mutex;
      [[maybe_unused]] std::lock_guard<std::mutex> guard(mutex);
      return LocalRaviartThomasFiniteElementFactory<D, d, R>::create(geometry_type, order);
    })
  {}
}; // ... LocalRaviartThomasFiniteElementFamily(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_RAVIART_THOMAS_HH
