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

#ifndef DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH
#define DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/raviart-thomas.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


/**
 * Applies
 * - scaling of shape functions to ensure basis*integrationElementNormal = 1
 * - flipping of shape functions to ensure continuity of normal component
 * - Piola transformation
 * - left-multiplication by the geometry transformations jacobian inverse transpose in jacobian
 */
template <class GL, class R = double>
class RaviartThomasGlobalBasis : public GlobalBasisInterface<GL, GL::dimension, 1, R>
{
  using ThisType = RaviartThomasGlobalBasis<GL, R>;
  using BaseType = GlobalBasisInterface<GL, GL::dimension, 1, R>;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::D;
  using typename BaseType::E;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalizedType;

  RaviartThomasGlobalBasis(const ThisType&) = default;
  RaviartThomasGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RaviartThomasGlobalBasis(const GridViewType& grid_view,
                           const int order,
                           const LocalRaviartThomasFiniteElementFamily<D, d, R>& local_finite_elements,
                           const FiniteVolumeMapper<GL>& element_indices,
                           const std::vector<DynamicVector<R>>& fe_data)
    : grid_view_(grid_view)
    , order_(order)
    , local_finite_elements_(local_finite_elements)
    , element_indices_(element_indices)
    , fe_data_(fe_data)
    , max_size_(0)
  {
    this->update_after_adapt();
  }

  size_t max_size() const override final
  {
    return max_size_;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedType> localize() const override final
  {
    return std::make_unique<LocalizedRaviartThomasGlobalBasis>(*this);
  }

  void update_after_adapt() override final
  {
    max_size_ = 0;
    for (const auto& gt : grid_view_.indexSet().types(0))
      max_size_ = std::max(max_size_, local_finite_elements_.get(gt, order_).size());
  }

private:
  class LocalizedRaviartThomasGlobalBasis : public LocalizedGlobalFiniteElementInterface<E, r, rC, R>
  {
    using ThisType = LocalizedRaviartThomasGlobalBasis;
    using BaseType = LocalizedGlobalFiniteElementInterface<E, r, rC, R>;

  public:
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::LocalFiniteElementType;
    using typename BaseType::RangeType;

    LocalizedRaviartThomasGlobalBasis(const RaviartThomasGlobalBasis<GL, R>& self)
      : BaseType()
      , self_(self)
      , set_data_in_post_bind_(true)
    {}

    LocalizedRaviartThomasGlobalBasis(const ThisType&) = default;
    LocalizedRaviartThomasGlobalBasis(ThisType&&) = default;

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
      current_local_fe_ = XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, d>>(
          self_.local_finite_elements_.get(elemnt.geometry().type(), self_.order_));
      if (set_data_in_post_bind_)
        current_fe_data_ = self_.fe_data_.at(self_.element_indices_.global_index(elemnt, 0));
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
      // evaluate shape functions
      current_local_fe_.access().basis().evaluate(point_in_reference_element, result);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < current_local_fe_.access().size(); ++ii)
        result[ii] *= current_fe_data_[ii];
      // apply piola transformation
      const auto J_T = this->element().geometry().jacobianTransposed(point_in_reference_element);
      const auto det_J_T = std::abs(J_T.determinant());
      RangeType tmp_value;
      for (size_t ii = 0; ii < current_local_fe_.access().size(); ++ii) {
        J_T.mtv(result[ii], tmp_value);
        tmp_value /= det_J_T;
        result[ii] = tmp_value;
      }
    }

    void jacobians(const DomainType& point_in_reference_element,
                   std::vector<DerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!current_local_fe_.valid(), Exceptions::not_bound_to_an_element_yet, "");
      this->assert_inside_reference_element(point_in_reference_element);
      // evaluate jacobian of shape functions
      current_local_fe_.access().basis().jacobian(point_in_reference_element, result);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < current_local_fe_.access().size(); ++ii)
        result[ii] *= current_fe_data_[ii];
      // apply piola transformation
      const auto J_T = this->element().geometry().jacobianTransposed(point_in_reference_element);
      const auto det_J_T = std::abs(J_T.determinant());
      const auto J_inv_T = this->element().geometry().jacobianInverseTransposed(point_in_reference_element);
      auto tmp_jacobian_row = result[0][0];
      for (size_t ii = 0; ii < current_local_fe_.access().size(); ++ii) {
        for (size_t jj = 0; jj < d; ++jj) {
          J_inv_T.mtv(result[ii][jj], tmp_jacobian_row);
          J_T.mtv(tmp_jacobian_row, result[ii][jj]);
          result[ii][jj] /= det_J_T;
        }
      }
    } // ... jacobian(...)

    // required by LocalizedGlobalFiniteElementInterface

    const LocalFiniteElementType& finite_element() const override final
    {
      return current_local_fe_.access();
    }

    DynamicVector<R> default_data(const GeometryType& geometry_type) const override final
    {
      return DynamicVector<R>(self_.local_finite_elements_.get(geometry_type, self_.order_).size(), 1.);
    }

    DynamicVector<R> backup() const override final
    {
      return current_fe_data_;
    }

    void restore(const ElementType& elemnt, const DynamicVector<R>& data) override final
    {
      set_data_in_post_bind_ = false;
      this->bind(elemnt);
      set_data_in_post_bind_ = true;
      DUNE_THROW_IF(data.size() != current_local_fe_.access().size(),
                    Exceptions::finite_element_error,
                    "data.size() = " << data.size()
                                     << "\ncurrent_local_fe_.access().size() = " << current_local_fe_.access().size());
      current_fe_data_ = data;
    } // ... restore(...)

    using BaseType::interpolate;

    void interpolate(const std::function<RangeType(DomainType)>& /*element_function*/,
                     const int /*order*/,
                     DynamicVector<R>& /*dofs*/) const override final
    {
      DUNE_THROW(NotImplemented, "");
    }

  private:
    const RaviartThomasGlobalBasis<GL, R>& self_;
    bool set_data_in_post_bind_;
    XT::Common::ConstStorageProvider<LocalFiniteElementInterface<D, d, R, d>> current_local_fe_;
    DynamicVector<R> current_fe_data_;
  }; // class LocalizedRaviartThomasGlobalBasis

  const GridViewType& grid_view_;
  const int order_;
  const LocalRaviartThomasFiniteElementFamily<D, d, R>& local_finite_elements_;
  const FiniteVolumeMapper<GL>& element_indices_;
  const std::vector<DynamicVector<R>>& fe_data_;
  size_t max_size_;
}; // class RaviartThomasGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH
