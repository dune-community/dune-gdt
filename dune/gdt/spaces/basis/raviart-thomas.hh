// This file is part of the dune-gdt project:
//   https://github.com/dune-comparamnity/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH
#define DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
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
  using typename BaseType::LocalizedBasisType;
  using typename BaseType::ShapeFunctionsType;
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC>;

  RaviartThomasGlobalBasis(const ThisType&) = default;
  RaviartThomasGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RaviartThomasGlobalBasis(
      const GridViewType& grd_vw,
      const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements,
      const std::shared_ptr<FiniteVolumeMapper<GL>> entity_indices,
      const std::shared_ptr<std::vector<std::vector<R>>> switches)
    : grid_view_(grd_vw)
    , finite_elements_(finite_elements)
    , element_indices_(entity_indices)
    , switches_(switches)
    , max_size_(0)
  {
    for (const auto& geometry_and_fe_pair : *finite_elements_)
      max_size_ = std::max(max_size_, geometry_and_fe_pair.second->basis().size());
  }

  virtual const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const override final
  {
    const auto finite_element_search_result = finite_elements_->find(geometry_type);
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   geometry_type = " << geometry_type);
    return finite_element_search_result->second->basis();
  }

  size_t max_size() const override final
  {
    return max_size_;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedBasisType> localize() const override final
  {
    return std::make_unique<LocalizedRaviartThomasGlobalBasis>(*this);
  }

private:
  class LocalizedRaviartThomasGlobalBasis : public XT::Functions::ElementFunctionSetInterface<E, r, rC, R>
  {
    using ThisType = LocalizedRaviartThomasGlobalBasis;
    using BaseType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

  public:
    using typename BaseType::DerivativeRangeType;
    using typename BaseType::DomainType;
    using typename BaseType::ElementType;
    using typename BaseType::RangeType;

    LocalizedRaviartThomasGlobalBasis(const RaviartThomasGlobalBasis<GL, R>& self)
      : BaseType()
      , self_(self)
      , shape_functions_(nullptr)
      , element_index_(0)
    {}

    LocalizedRaviartThomasGlobalBasis(const ThisType&) = default;
    LocalizedRaviartThomasGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return self_.max_size_;
    }

  protected:
    void post_bind(const ElementType& elemnt) override final
    {
      shape_functions_ = std::make_unique<XT::Common::ConstStorageProvider<ShapeFunctionsType>>(
          self_.shape_functions(elemnt.geometry().type()));
      element_index_ = self_.element_indices_->global_index(elemnt, 0);
      if ((*self_.switches_)[element_index_].size() != shape_functions_->access().size())
        DUNE_THROW(Exceptions::basis_error,
                   "shape_functions_->access().size() = " << shape_functions_->access().size()
                                                          << "\n   switches.size() = "
                                                          << (*self_.switches_)[element_index_].size());
    }

  public:
    size_t size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      return shape_functions_->access().size();
    }

    int order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      return shape_functions_->access().order();
    }

    using BaseType::evaluate;
    using BaseType::jacobians;

    void evaluate(const DomainType& point_in_reference_element,
                  std::vector<RangeType>& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      this->assert_inside_reference_element(point_in_reference_element);
      // evaluate shape functions
      shape_functions_->access().evaluate(point_in_reference_element, result);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < shape_functions_->access().size(); ++ii)
        result[ii] *= (*self_.switches_)[element_index_][ii];
      // apply piola transformation
      const auto J_T = this->element().geometry().jacobianTransposed(point_in_reference_element);
      const auto det_J_T = std::abs(J_T.determinant());
      RangeType tmp_value;
      for (size_t ii = 0; ii < shape_functions_->access().size(); ++ii) {
        J_T.mtv(result[ii], tmp_value);
        tmp_value /= det_J_T;
        result[ii] = tmp_value;
      }
    }

    void jacobians(const DomainType& point_in_reference_element,
                   std::vector<DerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      this->assert_inside_reference_element(point_in_reference_element);
      // evaluate jacobian of shape functions
      shape_functions_->access().jacobian(point_in_reference_element, result);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < shape_functions_->access().size(); ++ii)
        result[ii] *= (*self_.switches_)[element_index_][ii];
      // apply piola transformation
      const auto J_T = this->element().geometry().jacobianTransposed(point_in_reference_element);
      const auto det_J_T = std::abs(J_T.determinant());
      const auto J_inv_T = this->element().geometry().jacobianInverseTransposed(point_in_reference_element);
      auto tmp_jacobian_row = result[0][0];
      for (size_t ii = 0; ii < shape_functions_->access().size(); ++ii) {
        for (size_t jj = 0; jj < d; ++jj) {
          J_inv_T.mtv(result[ii][jj], tmp_jacobian_row);
          J_T.mtv(tmp_jacobian_row, result[ii][jj]);
          result[ii][jj] /= det_J_T;
        }
      }
    } // ... jacobian(...)

  private:
    const RaviartThomasGlobalBasis<GL, R>& self_;
    std::unique_ptr<XT::Common::ConstStorageProvider<ShapeFunctionsType>> shape_functions_;
    size_t element_index_;
  }; // class LocalizedRaviartThomasGlobalBasis

  const GridViewType& grid_view_;
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  const std::shared_ptr<FiniteVolumeMapper<GL>> element_indices_;
  const std::shared_ptr<std::vector<std::vector<R>>> switches_;
  size_t max_size_;
}; // class RaviartThomasGlobalBasis

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH
