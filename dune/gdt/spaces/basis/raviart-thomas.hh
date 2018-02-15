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

#include <dune/xt/grid/type_traits.hh>

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
  using typename BaseType::E;
  using typename BaseType::D;
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::ShapeFunctionsType;
  using typename BaseType::LocalizedBasisType;
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
    , entity_indices_(entity_indices)
    , switches_(switches)
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
    return std::make_unique<LocalizedRaviartThomasGlobalBasis>(
        element, shape_functions(element.geometry().type()), (*switches_)[entity_indices_->global_index(element, 0)]);
  }

private:
  class LocalizedRaviartThomasGlobalBasis : public XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>
  {
    using ThisType = LocalizedRaviartThomasGlobalBasis;
    using BaseType = XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>;

  public:
    using typename BaseType::EntityType;
    using typename BaseType::DomainType;
    using typename BaseType::RangeType;
    using typename BaseType::JacobianRangeType;

    LocalizedRaviartThomasGlobalBasis(const EntityType& elemnt,
                                      const ShapeFunctionsType& shape_funcs,
                                      const std::vector<R>& switches)
      : BaseType(elemnt)
      , shape_functions_(shape_funcs)
      , switches_(switches)
    {
      if (switches_.size() != shape_functions_.size())
        DUNE_THROW(
            Exceptions::basis_error,
            "shape_functions_.size() = " << shape_functions_.size() << "\n   switches.size() = " << switches_.size());
    }

    LocalizedRaviartThomasGlobalBasis(const ThisType&) = default;
    LocalizedRaviartThomasGlobalBasis(ThisType&&) = default;

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

    using BaseType::evaluate;

    void evaluate(const DomainType& xx,
                  std::vector<RangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      // evaluate shape functions
      ret = shape_functions_.evaluate(xx);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii)
        ret[ii] *= switches_[ii];
      // apply piola transformation
      const auto J_T = this->entity().geometry().jacobianTransposed(xx);
      const auto det_J_T = std::abs(J_T.determinant());
      RangeType tmp_value;
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii) {
        J_T.mtv(ret[ii], tmp_value);
        tmp_value /= det_J_T;
        ret[ii] = tmp_value;
      }
    }

    using BaseType::jacobian;

    void jacobian(const DomainType& xx,
                  std::vector<JacobianRangeType>& ret,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      // evaluate jacobian of shape functions
      ret = shape_functions_.jacobian(xx);
      // flip and scale shape functions to ensure
      // - continuity of normal component and
      // - basis*integrationElementNormal = 1
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii)
        ret[ii] *= switches_[ii];
      // apply piola transformation
      const auto J_T = this->entity().geometry().jacobianTransposed(xx);
      const auto det_J_T = std::abs(J_T.determinant());
      const auto J_inv_T = this->entity().geometry().jacobianInverseTransposed(xx);
      auto tmp_jacobian_row = ret[0][0];
      for (size_t ii = 0; ii < shape_functions_.size(); ++ii) {
        for (size_t jj = 0; jj < d; ++jj) {
          J_inv_T.mtv(ret[ii][jj], tmp_jacobian_row);
          J_T.mtv(tmp_jacobian_row, ret[ii][jj]);
          ret[ii][jj] /= det_J_T;
        }
      }
    } // ... jacobian(...)

  private:
    const ShapeFunctionsType& shape_functions_;
    const std::vector<R>& switches_;
  }; // class LocalizedRaviartThomasGlobalBasis

  const GridViewType& grid_view_;
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  const std::shared_ptr<FiniteVolumeMapper<GL>> entity_indices_;
  const std::shared_ptr<std::vector<std::vector<R>>> switches_;
}; // class RaviartThomasGlobalBasis

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_RAVIART_THOMAS_HH
