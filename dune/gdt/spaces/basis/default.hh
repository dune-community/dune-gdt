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

#ifndef DUNE_GDT_SPACES_BASIS_DEFAULT_HH
#define DUNE_GDT_SPACES_BASIS_DEFAULT_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


/**
 * Applies no transformation in evaluate, but left-multiplication by the geometry transformations jacobian inverse
 * transpose in jacobian.
 */
template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class DefaultGlobalBasis : public GlobalBasisInterface<GV, r, rC, R>
{
  using ThisType = DefaultGlobalBasis<GV, r, rC, R>;
  using BaseType = GlobalBasisInterface<GV, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::E;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalizedType;
  using FiniteElementFamilyType = LocalFiniteElementFamilyInterface<D, d, R, r, rC>;

  DefaultGlobalBasis(const ThisType&) = default;
  DefaultGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  DefaultGlobalBasis(const GridViewType& grid_view,
                     const FiniteElementFamilyType& local_finite_elements,
                     const int order)
    : grid_view_(grid_view)
    , local_finite_elements_(local_finite_elements)
    , order_(order)
    , max_size_(0)
  {}

  size_t max_size() const override final
  {
    return max_size_;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedType> localize() const override final
  {
    return std::make_unique<LocalizedDefaultGlobalBasis>(*this);
  }

  void update_after_adapt() override final
  {
    max_size_ = 0;
    for (const auto& gt : grid_view_.indexSet().types(0))
      max_size_ = std::max(max_size_, local_finite_elements_.get(gt, order_).size());
  }

private:
  class LocalizedDefaultGlobalBasis : public LocalizedGlobalFiniteElementInterface<E, r, rC, R>
  {
    using ThisType = LocalizedDefaultGlobalBasis;
    using BaseType = LocalizedGlobalFiniteElementInterface<E, r, rC, R>;
    static_assert(rC == 1,
                  "Not implemented for rC > 1, check that all functions (especially jacobians(...)) work properly "
                  "before removing this assert!");

  public:
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::LocalFiniteElementType;
    using typename BaseType::RangeType;

    LocalizedDefaultGlobalBasis(const DefaultGlobalBasis<GV, r, rC, R>& self)
      : BaseType()
      , self_(self)
    {}

    LocalizedDefaultGlobalBasis(const ThisType&) = delete;
    LocalizedDefaultGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    // required by XT::Functions::ElementFunctionSet

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return self_.max_size_;
    }

  protected:
    void post_bind(const ElementType& elemnt) override final
    {
      current_local_fe_ = XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, r, rC>>(
          self_.local_finite_elements_.get(elemnt.geometry().type(), self_.order_));
    }

  public:
    size_t size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!current_local_fe_.valid(), Exceptions::not_bound_to_an_element_yet, "");
      return current_local_fe_.access().size();
    }

    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!current_local_fe_.valid(), Exceptions::not_bound_to_an_element_yet, "");
      return current_local_fe_.access().basis().order();
    }

    using BaseType::evaluate;
    using BaseType::jacobians;

    void evaluate(const DomainType& point_in_reference_element,
                  std::vector<RangeType>& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!current_local_fe_.valid(), Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      current_local_fe_.access().basis().evaluate(point_in_reference_element, result);
    }

    void jacobians(const DomainType& point_in_reference_element,
                   std::vector<DerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!current_local_fe_.valid(), Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      // evaluate jacobian of shape functions
      current_local_fe_.access().basis().jacobian(point_in_reference_element, result);
      // Apply transformation:
      // Let f: E -> R^r be a basis function, and g: E' -> E be the mapping from reference to actual element, then f
      // \circ g is a shape function. We have the chain rule J_f = J(f \circ g \circ g^{-1}) = J(f \circ g) J_g^{-1}.
      // Applying the transpose to that equation gives
      // J_f^T = J_g^{-T} J(f \circ g)^T,
      // so we have to multiply J_inv_T from the left to the transposed shape function jacobians (i.e. the shape
      // function gradients) to get the transposed jacobian of the basis function (basis function gradient).
      const auto J_inv_T = this->element().geometry().jacobianInverseTransposed(point_in_reference_element);
      auto tmp_value = result[0][0];
      for (size_t ii = 0; ii < current_local_fe_.access().basis().size(); ++ii)
        for (size_t rr = 0; rr < r; ++rr) {
          J_inv_T.mv(result[ii][rr], tmp_value);
          result[ii][rr] = tmp_value;
        }
    } // ... jacobian(...)

    // required by LocalizedGlobalFiniteElementInterface

    const LocalFiniteElementType& finite_element() const override final
    {
      return current_local_fe_.access();
    }

    DynamicVector<R> default_data(const GeometryType& /*geometry_type*/) const override final
    {
      return DynamicVector<R>();
    }

    DynamicVector<R> backup() const override final
    {
      return DynamicVector<R>();
    }

    void restore(const ElementType& elemnt, const DynamicVector<R>& data) override final
    {
      this->bind(elemnt);
      DUNE_THROW_IF(data.size() != 0, Exceptions::finite_element_error, "data.size() = " << data.size());
    }

    using BaseType::interpolate;

    void interpolate(const std::function<RangeType(DomainType)>& element_function,
                     const int order,
                     DynamicVector<R>& dofs) const override final
    {
      current_local_fe_.access().interpolation().interpolate(element_function, order, dofs);
    }

  private:
    const DefaultGlobalBasis<GV, r, rC, R>& self_;
    XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, r, rC>> current_local_fe_;
  }; // class LocalizedDefaultGlobalBasis

  const GridViewType& grid_view_;
  const FiniteElementFamilyType& local_finite_elements_;
  const int order_;
  size_t max_size_;
}; // class DefaultGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_DEFAULT_HH
