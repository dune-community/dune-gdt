// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_SKELETON_HH
#define DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_SKELETON_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
#include <dune/gdt/local/finite-elements/lagrange.hh>

#include "interfaces.hh"
#include "finite-volume.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1>
class FiniteVolumeSkeletonMapper : public MapperInterface<GV>
{
  static_assert(r == 1, "The FiniteVolumeSkeletonMapper is not yet available for r > 1!");
  static_assert(rC == 1, "The FiniteVolumeSkeletonMapper is not yet available for rC > 1!");

  using ThisType = FiniteVolumeSkeletonMapper;
  using BaseType = MapperInterface<GV>;

public:
  using BaseType::d;
  using typename BaseType::D;

private:
  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV>;
  using LocalFiniteElementFamily = LocalLagrangeFiniteElementFamily<D, d, double, r>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  FiniteVolumeSkeletonMapper(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , local_finite_elements_()
    , element_indices_(grid_view_)
    , mapper_(grid_view_, mcmgLayout(Codim<1>()))
    , element_to_global_indices_of_intersections_map_(element_indices_.size())
    , max_local_size_(0)
  {
    this->update_after_adapt();
  }

  FiniteVolumeSkeletonMapper(const ThisType&) = default;
  FiniteVolumeSkeletonMapper(ThisType&&) = default;

  FiniteVolumeSkeletonMapper& operator=(const ThisType&) = delete;
  FiniteVolumeSkeletonMapper& operator=(ThisType&&) = delete;

  const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const LocalFiniteElementCoefficientsInterface<D, d>&
  local_coefficients(const GeometryType& geometry_type) const override final
  {
    DUNE_THROW(Exceptions::mapper_error, "Not implemented yet for a skeleton mapper!");
    return local_finite_elements_.get(geometry_type, 0).coefficients();
  }

  size_t size() const override final
  {
    return mapper_.size();
  }

  size_t max_local_size() const override final
  {
    return max_local_size_;
  }

  size_t local_size(const ElementType& element) const override final
  {
    return element.subEntities(1u);
  }

  size_t global_index(const ElementType& element, const size_t local_index) const override final
  {
    const auto& global_indices_of_intersections =
        element_to_global_indices_of_intersections_map_[element_indices_.global_index(element, 0)];
    DUNE_THROW_IF(local_index >= global_indices_of_intersections.size(),
                  Exceptions::mapper_error,
                  "local_size(element) = " << global_indices_of_intersections.size()
                                           << "\n   local_index = " << local_index);
    return global_indices_of_intersections[local_index];
  }

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    const auto& global_indices_of_intersections =
        element_to_global_indices_of_intersections_map_[element_indices_.global_index(element, 0)];
    const auto local_sz = global_indices_of_intersections.size();
    if (indices.size() < local_sz)
      indices.resize(local_sz);
    for (size_t ii = 0; ii < local_sz; ++ii)
      indices[ii] = global_indices_of_intersections[ii];
  } // ... global_indices(...)

  void update_after_adapt() override final
  {
    element_indices_.update_after_adapt();
    element_to_global_indices_of_intersections_map_.resize(element_indices_.size());
    mapper_.update();
    for (auto&& element : elements(grid_view_)) {
      const auto local_sz = this->local_size(element);
      max_local_size_ = std::max(max_local_size_, local_sz);
      auto& global_indices_of_intersections =
          element_to_global_indices_of_intersections_map_[element_indices_.global_index(element, 0)];
      global_indices_of_intersections.resize(local_sz);
      for (auto&& intersection : intersections(grid_view_, element)) {
        DUNE_THROW_IF(!intersection.conforming(),
                      Exceptions::mapper_error,
                      "Skeleton mapper not implemented for nonconforming intersections!");
        auto local_index = XT::Common::numeric_cast<size_t>(intersection.indexInInside());
        auto global_index = XT::Common::numeric_cast<size_t>(mapper_.index(element.template subEntity<1>(local_index)));
        global_indices_of_intersections[local_index] = global_index;
      }
    }
  } // ... update_after_adapt(...)

private:
  const GridViewType& grid_view_;
  LocalFiniteElementFamily local_finite_elements_;
  FiniteVolumeMapper<GV> element_indices_;
  Implementation mapper_;
  std::vector<DynamicVector<size_t>> element_to_global_indices_of_intersections_map_;
  size_t max_local_size_;
}; // class FiniteVolumeSkeletonMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_SKELETON_HH
