// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_SPACES_MAPPER_CONTINUOUS_HH
#define DUNE_GDT_SPACES_MAPPER_CONTINUOUS_HH

#include <set>

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/type_traits.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GV, class FiniteElement>
class ContinuousMapper : public MapperInterface<GV>
{
  static_assert(is_local_finite_element<FiniteElement>::value, "");
  using ThisType = ContinuousMapper<GV, FiniteElement>;
  using BaseType = MapperInterface<GV>;
  static const size_t r = FiniteElement::r;

  template <int d>
  struct GeometryTypeLayout
  {
    GeometryTypeLayout(std::set<GeometryType>&& types)
      : types_(std::move(types))
    {}

    GeometryTypeLayout(const GeometryTypeLayout<d>&) = default;
    GeometryTypeLayout(GeometryTypeLayout<d>&&) = default;

    bool contains(const GeometryType& gt) const
    {
      return types_.count(gt) > 0;
    }

    const std::set<GeometryType> types_;
  };

  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV, GeometryTypeLayout>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  ContinuousMapper(const GridViewType& grd_vw,
                   const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElement>>>& finite_elements)
    : grid_view_(grd_vw)
    , finite_elements_(finite_elements)
    , max_local_size_(0)
    , mapper_(nullptr)
  {
    DUNE_THROW_IF(d == 3 && finite_elements_->size() != 1, // Probably due to non-conforming intersections.
                  Exceptions::mapper_error,
                  "The mapper does not seem to work with multiple finite elements in 3d!");
    // collect all entities (for all codims) which are used to attach DoFs to
    for (auto&& geometry_type : grid_view_.indexSet().types(0)) {
      // get the finite element for this geometry type
      const auto finite_element_search_result = finite_elements_->find(geometry_type);
      DUNE_THROW_IF(finite_element_search_result == finite_elements_->end(),
                    Exceptions::mapper_error,
                    "Missing finite element for the required geometry type " << geometry_type << "!");
      const auto& finite_element = *finite_element_search_result->second;
      max_local_size_ = std::max(max_local_size_, finite_element.size());
      // loop over all keys of this finite element
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto& coeffs = finite_element.coefficients();
      for (size_t ii = 0; ii < coeffs.size(); ++ii) {
        const auto& local_key = coeffs.local_key(ii);
        // Currently only works if each subEntity has exactly r DoFs, if there is a variable number of DoFs per element
        // we would need to do more complicated things in the global index mapping.
        DUNE_THROW_IF(!(local_key.index() < r),
                      Exceptions::mapper_error,
                      "This case is not covered yet, when we have a variable number of DoFs per (sub)entity!");
        // find the (sub)entity for this key
        const auto sub_entity = local_key.subEntity();
        const auto codim = local_key.codim();
        const auto& subentity_geometry_type = reference_element.type(sub_entity, codim);
        // and add the respective geometry type
        all_DoF_attached_geometry_types_.insert(subentity_geometry_type);
      }
    }
    DUNE_THROW_IF(all_DoF_attached_geometry_types_.size() == 0,
                  Exceptions::mapper_error,
                  "This must not happen, the finite elements report no DoFs attached to (sub)entities!");
    mapper_ = std::make_shared<Implementation>(grid_view_, [&](const auto& geometry_type, const auto /*grid_dim*/) {
      return all_DoF_attached_geometry_types_.count(geometry_type) > 0;
    });
  } // ... ContinuousMapper(...)

  ContinuousMapper(const ThisType&) = default;
  ContinuousMapper(ThisType&&) = default;

  ContinuousMapper& operator=(const ThisType&) = delete;
  ContinuousMapper& operator=(ThisType&&) = delete;

  const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const LocalFiniteElementCoefficientsInterface<D, d>&
  local_coefficients(const GeometryType& geometry_type) const override final
  {
    const auto finite_element_search_result = finite_elements_->find(geometry_type);
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid view did not report all geometry types!"
                     << "\n   geometry_type = " << geometry_type);
    return finite_element_search_result->second->coefficients();
  }

  size_t size() const override final
  {
    return mapper_->size();
  }

  size_t max_local_size() const override final
  {
    return max_local_size_;
  }

  size_t local_size(const ElementType& element) const override final
  {
    return local_coefficients(element.geometry().type()).size();
  }

  size_t global_index(const ElementType& element, const size_t local_index) const override final
  {
    const auto& coeffs = local_coefficients(element.geometry().type());
    if (local_index >= coeffs.size())
      DUNE_THROW(Exceptions::mapper_error,
                 "local_size(element) = " << coeffs.size() << "\n   local_index = " << local_index);
    const auto& local_key = coeffs.local_key(local_index);
    return r * mapper_->subIndex(element, local_key.subEntity(), local_key.codim()) + local_key.index();
  } // ... mapToGlobal(...)

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    const auto& coeffs = local_coefficients(element.geometry().type());
    const auto local_sz = coeffs.size();
    if (indices.size() < local_sz)
      indices.resize(local_sz, 0);
    for (size_t ii = 0; ii < local_sz; ++ii) {
      const auto& local_key = coeffs.local_key(ii);
      // No need to assert local_key.index() == 0, has been checked in the ctor!
      indices[ii] = r * mapper_->subIndex(element, local_key.subEntity(), local_key.codim()) + local_key.index();
    }
  } // ... globalIndices(...)

private:
  const GridViewType& grid_view_;
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElement>>> finite_elements_;
  size_t max_local_size_;
  std::set<GeometryType> all_DoF_attached_geometry_types_;
  std::shared_ptr<Implementation> mapper_;
}; // class ContinuousMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_CONTINUOUS_HH
