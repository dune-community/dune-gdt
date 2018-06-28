// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)

#ifndef DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_HH
#define DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_HH

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
#include <dune/gdt/local/finite-elements/lagrange.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1>
class FiniteVolumeMapper : public MapperInterface<GV>
{
  static_assert(rC == 1, "The FiniteVolumeMapper is not yet available for rC > 1!");

  using ThisType = FiniteVolumeMapper<GV, r, rC>;
  using BaseType = MapperInterface<GV>;

public:
  using typename BaseType::D;
  using BaseType::d;

private:
  template <int>
  struct Codim0EntityFilter
  {
    bool contains(const GeometryType& gt) const
    {
      return gt.dim() == d;
    }
  };

  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV, Codim0EntityFilter>;
  using FiniteElement = LocalFiniteElementInterface<D, d, double, r, rC>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  FiniteVolumeMapper(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , mapper_(new Implementation(grid_view_))
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElement>>())
  {
    // create finite elements (we only need the coefficients)
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_local_lagrange_finite_element<D, d, double, r>(geometry_type, 0)));
  }

  FiniteVolumeMapper(const ThisType&) = default;
  FiniteVolumeMapper(ThisType&&) = default;

  FiniteVolumeMapper& operator=(const ThisType&) = delete;
  FiniteVolumeMapper& operator=(ThisType&&) = delete;

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
                     << "\n   geometry_type = "
                     << geometry_type);
    return finite_element_search_result->second->coefficients();
  } // ... local_coefficients(...)

  size_t size() const override final
  {
    return mapper_->size() * r * rC;
  }

  size_t max_local_size() const override final
  {
    return r * rC;
  }

  size_t local_size(const ElementType& /*element*/) const override final
  {
    return r * rC;
  }

  size_t global_index(const ElementType& element, const size_t local_index) const override final
  {
    if (local_index >= r * rC)
      DUNE_THROW(Exceptions::mapper_error, "local_size(element) = " << r * rC << "\n   local_index = " << local_index);
    return mapper_->index(element) * r * rC + local_index;
  }

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    if (indices.size() < r * rC)
      indices.resize(r * rC);
    size_t local_index = 0;
    for (size_t ii = 0; ii < r; ++ii)
      for (size_t jj = 0; jj < rC; ++jj) {
        indices[local_index] = mapper_->index(element) * r * rC + local_index;
        ++local_index;
      }
  } // ... global_indices(...)

private:
  const GridViewType& grid_view_;
  const std::shared_ptr<Implementation> mapper_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElement>>> finite_elements_;
}; // class FiniteVolumeMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_HH
