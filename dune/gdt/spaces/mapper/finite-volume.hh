// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

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

  using ThisType = FiniteVolumeMapper;
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

  FiniteVolumeMapper(const GridViewType& grd_vw)
    : grid_view_(grd_vw)
    , local_finite_elements_()
    , mapper_(grid_view_, mcmgElementLayout())
  {
    this->update_after_adapt();
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
    return local_finite_elements_.get(geometry_type, 0).coefficients();
  }

  size_t size() const override final
  {
    return mapper_.size() * r * rC;
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
    return mapper_.index(element) * r * rC + local_index;
  }

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    constexpr size_t local_sz = r * rC;
    if (indices.size() < local_sz)
      indices.resize(local_sz);
    size_t* indices_ptr = &indices[0];
    std::iota(indices_ptr, indices_ptr + local_sz, mapper_.index(element) * local_sz);
  } // ... global_indices(...)

  void update_after_adapt() override final
  {
    mapper_.update();
  }

private:
  const GridViewType& grid_view_;
  LocalFiniteElementFamily local_finite_elements_;
  Implementation mapper_;
}; // class FiniteVolumeMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_FINITE_VOLUME_HH
