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
  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  DiscontinuousMapper(const GridViewType& grd_vw,
                      const LocalFiniteElementFamily& local_finite_elements,
                      const int order)
    : grid_view_(grd_vw)
    , local_finite_elements_(local_finite_elements)
    , order_(order)
    , mapper_(grid_view_, mcmgElementLayout())
    , offset_(mapper_.size())
  {
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
    return size_;
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
    if (local_index >= local_size(element))
      DUNE_THROW(Exceptions::mapper_error,
                 "local_size(element) = " << local_size(element) << "\n   local_index = " << local_index);
    return offset_[mapper_.index(element)] + local_index;
  }

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    const size_t offset = offset_[mapper_.index(element)];
    const auto local_sz = local_size(element);
    if (indices.size() < local_sz)
      indices.resize(local_sz, 0);
    for (size_t ii = 0; ii < local_sz; ++ii)
      indices[ii] = offset + ii;
  } // ... global_indices(...)

  void update_after_adapt() override final
  {
    mapper_.update();
    size_ = 0;
    max_num_dofs_ = 0;
    if (offset_.size() != mapper_.size())
      offset_.resize(mapper_.size());
    for (auto&& element : elements(grid_view_)) {
      offset_[mapper_.index(element)] = size_;
      const auto local_sz = this->local_size(element);
      size_ += local_sz;
      max_num_dofs_ = std::max(max_num_dofs_, local_sz);
    }
  }

private:
  const GridViewType& grid_view_;
  const LocalFiniteElementFamily& local_finite_elements_;
  const int order_;
  Implementation mapper_;
  std::vector<size_t> offset_;
  size_t max_num_dofs_;
  size_t size_;
}; // class DiscontinuousMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH
