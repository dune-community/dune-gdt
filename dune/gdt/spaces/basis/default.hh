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

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/functions/interfaces/local-functions.hh>

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
  {
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

  std::unique_ptr<LocalizedBasisType> localize(const ElementType& element) const override final
  {
    return std::make_unique<LocalizedDefaultGlobalBasis>(element, shape_functions(element.geometry().type()));
  }

private:
  class LocalizedDefaultGlobalBasis : public XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>
  {
    using ThisType = LocalizedDefaultGlobalBasis;
    using BaseType = XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>;

  public:
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    LocalizedDefaultGlobalBasis(const ElementType& elemnt, const ShapeFunctionsType& shape_funcs)
      : BaseType(elemnt)
      , shape_functions_(shape_funcs)
    {
    }

    LocalizedDefaultGlobalBasis(const ThisType&) = default;
    LocalizedDefaultGlobalBasis(ThisType&&) = default;

    ThisType& operator=(const ThisType&) = delete;
    ThisType& operator=(ThisType&&) = delete;

    size_t size() const override final
    {
      return shape_functions_.size();
    }

    size_t order(const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return XT::Common::numeric_cast<size_t>(shape_functions_.order());
    }

    void evaluate(const DomainType& xx,
                  std::vector<RangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret = shape_functions_.evaluate(xx);
    }

    std::vector<RangeType> evaluate(const DomainType& xx,
                                    const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      return shape_functions_.evaluate(xx);
    }

    void jacobian(const DomainType& xx,
                  std::vector<JacobianRangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      // evaluate jacobian of shape functions
      ret = shape_functions_.jacobian(xx);
      // apply transformation
      const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
      auto tmp_value = ret[0][0];
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii) {
        J_inv_T.mv(ret[ii][0], tmp_value);
        ret[ii][0] = tmp_value;
      }
    } // ... jacobian(...)

    std::vector<JacobianRangeType> jacobian(const DomainType& xx,
                                            const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      // evaluate jacobian of shape functions
      auto ret = shape_functions_.jacobian(xx);
      // apply transformation
      const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
      auto tmp_value = ret[0][0];
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii) {
        J_inv_T.mv(ret[ii][0], tmp_value);
        ret[ii][0] = tmp_value;
      }
      return ret;
    } // ... jacobian(...)

  private:
    const ShapeFunctionsType& shape_functions_;
  }; // class LocalizedDefaultGlobalBasis

  const GridViewType& grid_view_;
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
}; // class DefaultGlobalBasis


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_DEFAULT_HH
