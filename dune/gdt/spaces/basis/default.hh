// This file is part of the dune-gdt project:
//   https://github.com/dune-comparamnity/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_BASIS_DEFAULT_HH
#define DUNE_GDT_SPACES_BASIS_DEFAULT_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/interfaces/local-functions.hh>

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
  using typename BaseType::E;
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::ShapeFunctionsType;
  using typename BaseType::LocalizedBasisType;
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC>;

  DefaultGlobalBasis(const ThisType&) = default;
  DefaultGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  DefaultGlobalBasis(const GridViewType& grd_vw,
                     const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements)
    : grid_view_(grd_vw)
    , finite_elements_(finite_elements)
    , max_size_()
  {
    for (const auto& geometry_and_fe_pair : *finite_elements_)
      max_size_ = std::max(max_size_, geometry_and_fe_pair.second->basis().size());
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const override final
  {
    const auto finite_element_search_result = finite_elements_->find(geometry_type);
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   geometry_type = "
                     << geometry_type);
    return finite_element_search_result->second->basis();
  }

  std::unique_ptr<LocalizedBasisType> localize() const override final
  {
    return std::make_unique<LocalizedDefaultGlobalBasis>(*this);
  }

  std::unique_ptr<LocalizedBasisType> localize(const ElementType& element) const override final
  {
    return std::make_unique<LocalizedDefaultGlobalBasis>(*this, element);
  }

private:
  class LocalizedDefaultGlobalBasis : public XT::Functions::LocalFunctionSetInterface<E, r, rC, R>
  {
    using ThisType = LocalizedDefaultGlobalBasis;
    using BaseType = XT::Functions::LocalFunctionSetInterface<E, r, rC, R>;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::DerivativeRangeType;

    LocalizedDefaultGlobalBasis(const DefaultGlobalBasis<GV, r, rC, R>& self)
      : BaseType()
      , self_(self)
      , shape_functions_(nullptr)
    {
    }

    LocalizedDefaultGlobalBasis(const DefaultGlobalBasis<GV, r, rC, R>& self, const EntityType& elemnt)
      : BaseType(elemnt)
      , self_(self)
      , shape_functions_(nullptr)
    {
      post_bind(elemnt);
    }

    LocalizedDefaultGlobalBasis(const ThisType&) = default;
    LocalizedDefaultGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t max_size(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return self_.max_size_;
    }

  protected:
    void post_bind(const EntityType& elemnt) override final
    {
      shape_functions_ = std::make_unique<XT::Common::ConstStorageProvider<ShapeFunctionsType>>(
          self_.shape_functions(elemnt.geometry().type()));
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

    void evaluate(const DomainType& point_in_reference_element,
                  std::vector<RangeType>& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      this->assert_inside_reference_element(point_in_reference_element);
      shape_functions_->access().evaluate(point_in_reference_element, result);
    }

    void jacobians(const DomainType& point_in_reference_element,
                   std::vector<DerivativeRangeType>& result,
                   const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      DUNE_THROW_IF(!shape_functions_, Exceptions::not_bound_to_an_element_yet, "you need to call bind() first!");
      this->assert_inside_reference_element(point_in_reference_element);
      // evaluate jacobian of shape functions
      shape_functions_->access().jacobian(point_in_reference_element, result);
      // apply transformation
      const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(point_in_reference_element);
      auto tmp_value = result[0][0];
      for (size_t ii = 0; ii < shape_functions_->access().size(); ++ii) {
        J_inv_T.mv(result[ii][0], tmp_value);
        result[ii][0] = tmp_value;
      }
    } // ... jacobian(...)

  private:
    const DefaultGlobalBasis<GV, r, rC, R>& self_;
    std::unique_ptr<XT::Common::ConstStorageProvider<ShapeFunctionsType>> shape_functions_;
  }; // class LocalizedDefaultGlobalBasis

  const GridViewType& grid_view_;
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  size_t max_size_;
}; // class DefaultGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_DEFAULT_HH
