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
 *
 * \todo Somehow make use of the LocalRaviartThomasInterpolation and get rid of the duplicate code.
 */
template <class GL, class R = double>
class RaviartThomasGlobalBasis : public GlobalBasisInterface<GL, GL::dimension, 1, R>
{
  using ThisType = RaviartThomasGlobalBasis;
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

  RaviartThomasGlobalBasis(const GridViewType& grid_view,
                           const int order,
                           const LocalRaviartThomasFiniteElementFamily<D, d, R>& local_finite_elements,
                           const FiniteVolumeMapper<GL>& element_indices,
                           const std::vector<DynamicVector<R>>& fe_data,
                           const std::string& logging_prefix)
    : BaseType(logging_prefix.empty() ? "gdt" : "gdt.spaces.rt",
               logging_prefix.empty() ? "RtGlobalBasis" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , grid_view_(grid_view)
    , order_(order)
    , local_finite_elements_(local_finite_elements)
    , element_indices_(element_indices)
    , fe_data_(fe_data)
    , max_size_(0)
  {
    LOG_(info) << this->logging_id << "(grid_view=" << &grid_view << ", order=" << order
               << ", local_finite_elements=" << &local_finite_elements
               << ",\n   element_indices.size()=" << element_indices.size() << ", fe_data.size()=" << fe_data.size()
               << ")" << std::endl;
    this->update_after_adapt();
  }

  RaviartThomasGlobalBasis(const ThisType&) = default;

  RaviartThomasGlobalBasis(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;

  ThisType& operator=(ThisType&&) = delete;

  size_t max_size() const override final
  {
    return max_size_;
  }

  using BaseType::localize;

  std::unique_ptr<LocalizedType> localize() const override final
  {
    std::string derived_logging_prefix = "";
    if (this->logger.info_enabled) {
      derived_logging_prefix = this->logging_id + "_localized";
      this->logger.info() << this->logging_id << ".localize()" << std::endl;
    }
    return std::make_unique<LocalizedRaviartThomasGlobalBasis>(*this, derived_logging_prefix);
  }

  void update_after_adapt() override final
  {
    LOG_(info) << this->logging_id << ".update_after_adapt()" << std::endl;
    max_size_ = 0;
    for (const auto& gt : grid_view_.indexSet().types(0))
      max_size_ = std::max(max_size_, local_finite_elements_.get(gt, order_).size());
    LOG_(info) << "  max_size_ = " << max_size_ << std::endl;
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

    LocalizedRaviartThomasGlobalBasis(const RaviartThomasGlobalBasis<GL, R>& self,
                                      const std::string& logging_prefix = "")
      : BaseType(logging_prefix.empty() ? "gdt" : "gdt.spaces.rt",
                 logging_prefix.empty() ? "LocalizedRtGlobalBasis" : logging_prefix,
                 /*logging_disabled=*/logging_prefix.empty())
      , self_(self)
      , set_data_in_post_bind_(true)
    {
      LOG_(info) << this->logging_id << "(self=" << &self
                 << ", self.element_indices_.size()=" << self.element_indices_.size()
                 << ", self.fe_data_.size()=" << self.fe_data_.size() << ")" << std::endl;
    }

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
          self_.local_finite_elements_.get(elemnt.type(), self_.order_));
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

    void interpolate(const std::function<RangeType(DomainType)>& element_function,
                     const int element_function_order,
                     DynamicVector<R>& dofs) const override final
    {
      // prepare
      const size_t sz = this->size();
      if (dofs.size() != sz)
        dofs.resize(sz);
      std::vector<bool> local_key_was_handled(sz, false);
      // determine the face dofs, therefore walk the intersections
      for (auto&& intersection : intersections(self_.grid_view_, this->element())) {
        const auto intersection_to_local_key_map = finite_element().coefficients().local_key_indices(1);
        const auto intersection_index = intersection.indexInInside();
        const auto& local_keys_assosiated_with_intersection = intersection_to_local_key_map[intersection_index];
        if (local_keys_assosiated_with_intersection.size() > 0) {
          const auto intersection_fe =
              make_local_orthonormal_finite_element<D, d - 1, R>(intersection.type(), finite_element().order());
          const auto& intersection_Pk_basis = intersection_fe->basis();
          DUNE_THROW_IF(intersection_Pk_basis.size() != local_keys_assosiated_with_intersection.size(),
                        Exceptions::interpolation_error,
                        "intersection_Pk_basis.size() = " << intersection_Pk_basis.size()
                                                          << "\n   local_keys_assosiated_with_intersection.size() = "
                                                          << local_keys_assosiated_with_intersection.size());
          XT::LA::CommonDenseMatrix<R> lhs(
              local_keys_assosiated_with_intersection.size(), intersection_Pk_basis.size(), 0);
          XT::LA::CommonDenseVector<R> rhs(intersection_Pk_basis.size(), 0);
          // do a face quadrature
          for (auto&& quadrature_point :
               QuadratureRules<D, d - 1>::rule(intersection.type(),
                                               std::max(this->order() + intersection_Pk_basis.order(),
                                                        element_function_order + intersection_Pk_basis.order()))) {
            const auto point_on_reference_intersection = quadrature_point.position();
            const auto point_in_reference_element =
                intersection.geometryInInside().global(point_on_reference_intersection);
            const auto quadrature_weight = quadrature_point.weight();
            const auto normal = intersection.unitOuterNormal(point_on_reference_intersection);
            const auto integration_factor = intersection.geometry().integrationElement(point_on_reference_intersection);
            const auto rt_basis_values = this->evaluate_set(point_in_reference_element);
            const auto intersection_Pk_basis_values = intersection_Pk_basis.evaluate(point_on_reference_intersection);
            const auto element_function_values = element_function(point_in_reference_element);
            for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
              const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
              for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                lhs.add_to_entry(ii,
                                 jj,
                                 quadrature_weight * integration_factor * (rt_basis_values[local_key_index] * normal)
                                     * intersection_Pk_basis_values[jj]);
            }
            for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
              rhs[jj] += quadrature_weight * integration_factor * (element_function_values * normal)
                         * intersection_Pk_basis_values[jj];
          }
          XT::LA::CommonDenseVector<R> intersection_dofs(local_keys_assosiated_with_intersection.size(), 0);
          try {
            intersection_dofs = XT::LA::solve(lhs, rhs);
          } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
            DUNE_THROW(Exceptions::interpolation_error,
                       "Failed to solve for DoFs associated with intersection "
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
      const auto element_to_local_key_map = finite_element().coefficients().local_key_indices(0);
      const auto& local_keys_assosiated_with_element = element_to_local_key_map[0];
      if (local_keys_assosiated_with_element.size() > 0) {
        DUNE_THROW_IF(this->order() < 1,
                      Exceptions::interpolation_error,
                      "DoFs associated with the element only make sense for orders >= 1!");
        const auto element_fe =
            make_local_orthonormal_finite_element<D, d, R, d>(this->element().type(), finite_element().order() - 1);
        const auto& element_Pkminus1_basis = element_fe->basis();
        DUNE_THROW_IF(element_Pkminus1_basis.size() != local_keys_assosiated_with_element.size(),
                      Exceptions::interpolation_error,
                      "element_Pkminus1_basis.size() = " << element_Pkminus1_basis.size()
                                                         << "\n   local_keys_assosiated_with_element.size() = "
                                                         << local_keys_assosiated_with_element.size());
        XT::LA::CommonDenseMatrix<R> lhs(local_keys_assosiated_with_element.size(), element_Pkminus1_basis.size(), 0);
        XT::LA::CommonDenseVector<R> rhs(element_Pkminus1_basis.size(), 0);
        // do a volume quadrature
        for (auto&& quadrature_point :
             QuadratureRules<D, d>::rule(this->element().type(),
                                         std::max(this->order() + element_Pkminus1_basis.order(),
                                                  element_function_order + element_Pkminus1_basis.order()))) {
          const auto point_in_reference_element = quadrature_point.position();
          const auto quadrature_weight = quadrature_point.weight();
          const auto integration_factor = this->element().geometry().integrationElement(point_in_reference_element);
          const auto rt_basis_values = this->evaluate_set(point_in_reference_element);
          const auto element_Pkminus1_basis_values = element_Pkminus1_basis.evaluate(point_in_reference_element);
          const auto element_function_values = element_function(point_in_reference_element);
          for (size_t ii = 0; ii < local_keys_assosiated_with_element.size(); ++ii) {
            const size_t local_key_index = local_keys_assosiated_with_element[ii];
            for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
              lhs.add_to_entry(ii,
                               jj,
                               quadrature_weight * integration_factor
                                   * (rt_basis_values[local_key_index] * element_Pkminus1_basis_values[jj]));
          }
          for (size_t jj = 0; jj < element_Pkminus1_basis.size(); ++jj)
            rhs[jj] +=
                quadrature_weight * integration_factor * (element_function_values * element_Pkminus1_basis_values[jj]);
        }
        XT::LA::CommonDenseVector<R> element_dofs(local_keys_assosiated_with_element.size(), 0);
        try {
          element_dofs = XT::LA::solve(lhs, rhs);
        } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
          DUNE_THROW(Exceptions::interpolation_error,
                     "Failed to solve for volume DoFs, this was the original error:\n   " << ee.what());
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
        DUNE_THROW_IF(!local_key_was_handled[ii],
                      Exceptions::interpolation_error,
                      "The following DoF is neither associated with an intersection, nor with the element!"
                          << "\n   local DoF index: " << ii
                          << "\n   associated local_key: " << finite_element().coefficients().local_key(ii));
    } // ... interpolate(...)

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
