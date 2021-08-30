// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH
#define DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/lagrange.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GV, class LocalFiniteElementFamily>
class DiscontinuousMapper : public MapperInterface<GV>
{
  static_assert(is_local_finite_element_family<LocalFiniteElementFamily>::value, "");
  using ThisType = DiscontinuousMapper;
  using BaseType = MapperInterface<GV>;

public:
  using BaseType::d;
  using typename BaseType::D;

private:
  static constexpr size_t r = LocalFiniteElementFamily::r;
  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  DiscontinuousMapper(const GridViewType& grd_vw,
                      const LocalFiniteElementFamily& local_finite_elements,
                      const int order,
                      const bool dimwise_global_numbering = false)
    : grid_view_(grd_vw)
    , local_finite_elements_(local_finite_elements)
    , order_(order)
    , dimwise_global_numbering_((r == 1) ? false : dimwise_global_numbering) // does not make sense in the scalar case
    , mapper_(grid_view_, mcmgElementLayout())
    , offset_()
  {
    if (!dimwise_global_numbering_) // only use first component
      offset_[0].resize(mapper_.size());
    else
      for (size_t s = 0; s < r; ++s)
        offset_[s].resize(mapper_.size());
    this->update_after_adapt();
  }

  DiscontinuousMapper(const ThisType&) = default;
  DiscontinuousMapper(ThisType&&) = default;

  DiscontinuousMapper& operator=(const ThisType&) = delete;
  DiscontinuousMapper& operator=(ThisType&&) = delete;

  const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const LocalFiniteElementCoefficientsInterface<D, d>&
  local_coefficients(const GeometryType& geometry_type) const override final
  {
    return local_finite_elements_.get(geometry_type, order_).coefficients();
  }

  size_t size() const override final
  {
    return std::accumulate(size_.begin(), size_.end(), 0);
  }

  size_t max_local_size() const override final
  {
    return max_num_dofs_;
  }

  size_t local_size(const ElementType& element) const override final
  {
    return local_coefficients(element.type()).size();
  }

  size_t global_index(const ElementType& element, const size_t local_index) const override final
  {
    const auto local_sz = local_size(element);
    DUNE_THROW_IF(local_index >= local_sz,
                  Exceptions::mapper_error,
                  "local_size(element) = " << local_sz << "\n   local_index = " << local_index);
    if (!dimwise_global_numbering_)
      return offset_[0][mapper_.index(element)] + local_index;
    else {
      assert(local_sz % r == 0 && "This should not happen, see update_after_adapt()!");
      const size_t num_local_indices_per_range_dim = local_sz / r;
      const size_t dim_range_associated_with_local_index = int(local_index / num_local_indices_per_range_dim);
      const size_t local_index_within_dim_range = local_index % num_local_indices_per_range_dim;
      size_t base_offset = 0;
      for (size_t s = 0; s < dim_range_associated_with_local_index; ++s)
        base_offset += size_[s];
      return base_offset + offset_[dim_range_associated_with_local_index][mapper_.index(element)]
             + local_index_within_dim_range;
    }
  } // ... global_index(...)

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    const auto local_sz = local_size(element);
    if (indices.size() < local_sz)
      indices.resize(local_sz, 0);
    if (!dimwise_global_numbering_) {
      const size_t offset = offset_[0][mapper_.index(element)];
      for (size_t ii = 0; ii < local_sz; ++ii)
        indices[ii] = offset + ii;
    } else {
      assert(local_sz % r == 0 && "This should not happen, see update_after_adapt()!");
      const size_t num_local_indices_per_range_dim = local_sz / r;
      for (size_t local_index = 0; local_index < local_sz; ++local_index) {
        const size_t dim_range_associated_with_local_index = int(local_index / num_local_indices_per_range_dim);
        const size_t local_index_within_dim_range = local_index % num_local_indices_per_range_dim;
        size_t offset = 0;
        for (size_t s = 0; s < dim_range_associated_with_local_index; ++s)
          offset += size_[s];
        offset += offset_[dim_range_associated_with_local_index][mapper_.index(element)];
        indices[local_index] = offset + local_index_within_dim_range;
      }
    }
  } // ... global_indices(...)

  void update_after_adapt() override final
  {
    mapper_.update();
    max_num_dofs_ = 0;
    std::fill(size_.begin(), size_.end(), 0);
    auto prepare = [&](const auto& i) {
      if (offset_[i].size() != mapper_.size())
        offset_[i].resize(mapper_.size());
    };
    if (!dimwise_global_numbering_) // only use first component
      prepare(0);
    else
      for (size_t s = 0; s < r; ++s)
        prepare(s);
    for (auto&& element : elements(grid_view_)) {
      const auto local_sz = this->local_size(element);
      if (!dimwise_global_numbering_) {
        offset_[0][mapper_.index(element)] = size_[0];
        size_[0] += local_sz;
      } else {
        DUNE_THROW_IF(!local_finite_elements_.get(element.geometry().type(), order_).powered(),
                      Exceptions::mapper_error,
                      "Using dimwise numbering only makes sense for powered FEs!");
        DUNE_THROW_IF(local_sz % r != 0,
                      Exceptions::mapper_error,
                      "Not a power FE, as number of local keys associated with an element not a multiple of "
                      "dim_range:\n   this->local_size(element) = "
                          << local_sz << "\n   "
                          << "dim_range = " << size_t(r));
        for (size_t s = 0; s < r; ++s) {
          offset_[s][mapper_.index(element)] = size_[s];
          size_[s] += local_sz / r;
        }
      }
      max_num_dofs_ = std::max(max_num_dofs_, local_sz);
    }
  } // ... update_after_adapt(...)

private:
  const GridViewType& grid_view_;
  const LocalFiniteElementFamily& local_finite_elements_;
  const int order_;
  const bool dimwise_global_numbering_;
  Implementation mapper_;
  std::array<std::vector<size_t>, r> offset_;
  std::array<size_t, r> size_;
  size_t max_num_dofs_;
}; // class DiscontinuousMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH
