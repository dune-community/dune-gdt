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


template <class GV, class LocalFiniteElementFamily>
class ContinuousMapper : public MapperInterface<GV>
{
  static_assert(is_local_finite_element_family<LocalFiniteElementFamily>::value, "");
  using ThisType = ContinuousMapper;
  using BaseType = MapperInterface<GV>;

  using Implementation = MultipleCodimMultipleGeomTypeMapper<GV>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;

  ContinuousMapper(const GridViewType& grd_vw,
                   const LocalFiniteElementFamily& local_finite_elements,
                   const int fe_order)
    : grid_view_(grd_vw)
    , local_finite_elements_(local_finite_elements)
    , fe_order_(fe_order)
    , max_local_size_(0)
    , mapper_(grid_view_, [&](const auto& geometry_type, const auto /*grid_dim*/) {
      return (all_DoF_attached_geometry_types_.count(geometry_type) > 0) ? geometry_type_to_local_size_[geometry_type]
                                                                         : 0;
    })
  {
    if (d >= 2 && fe_order_ >= 3 && !XT::Grid::is_cube_alugrid<typename GV::Grid>::value
        && !XT::Grid::is_yaspgrid<typename GV::Grid>::value && !XT::Grid::is_uggrid<typename GV::Grid>::value)
      DUNE_THROW(Dune::NotImplemented,
                 "For order > 2, there are problems with the local-to-global mapping on some grids, see the comment in "
                 "the global_index method!");
    this->update_after_adapt();
  }

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
    return local_finite_elements_.get(geometry_type, fe_order_).coefficients();
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
    return local_coefficients(element.type()).size();
  }

  size_t global_index(const ElementType& element, const size_t local_index) const override final
  {
    const auto& coeffs = local_coefficients(element.type());
    if (local_index >= coeffs.size())
      DUNE_THROW(Exceptions::mapper_error,
                 "local_size(element) = " << coeffs.size() << "\n   local_index = " << local_index);
    const auto& local_key = coeffs.local_key(local_index);
    // TODO: If there are several DoFs on one subEntity (which is the case e.g. for third order lagrange elements), this
    // mapping only works if the DoFs are numbered consistently between the elements. For example, if there are two DoFs
    // on a shared edge between to codim 0 elements, the local_key.index() has to be the same for the same DoF in both
    // elements. This does not seem to be the case for the simplex grids, the same DoF might be assigned the index 0 in
    // one element and index 1 in the other element. Fixing this could be done by assigning an orientation to the edge
    // by looking at the (indices of the) vertices of the edge and reordering the local indices if the orientation is
    // not the same in all elements sharing the subentity.
    DEBUG_THROW_IF(d >= 2 && fe_order_ >= 3 && element.type() != Dune::GeometryTypes::cube(d),
                   Dune::NotImplemented,
                   "Not implemented for this element, see comment above!");
    return mapper_.subIndex(element, local_key.subEntity(), local_key.codim()) + local_key.index();
  } // ... mapToGlobal(...)

  using BaseType::global_indices;

  void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const override final
  {
    const auto& coeffs = local_coefficients(element.type());
    const auto local_sz = coeffs.size();
    assert(local_sz <= max_local_size_);
    if (indices.size() < local_sz)
      indices.resize(local_sz, 0);
    for (size_t ii = 0; ii < local_sz; ++ii) {
      const auto& local_key = coeffs.local_key(ii);
      indices[ii] = mapper_.subIndex(element, local_key.subEntity(), local_key.codim()) + local_key.index();
    }
  } // ... globalIndices(...)

  void update_after_adapt() override final
  {
    // Probably due to non-conforming intersections.
    DUNE_THROW_IF(d == 3 && grid_view_.indexSet().types(0).size() != 1,
                  Exceptions::mapper_error,
                  "The mapper does not seem to work with multiple finite elements in 3d!");
    // collect all entities (for all codims) which are used to attach DoFs to
    all_DoF_attached_geometry_types_.clear();
    for (auto&& geometry_type : grid_view_.indexSet().types(0)) {
      const auto& finite_element = local_finite_elements_.get(geometry_type, fe_order_);
      max_local_size_ = std::max(max_local_size_, finite_element.size());
      // loop over all keys of this finite element
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto& coeffs = finite_element.coefficients();
      const auto& local_key_indices = coeffs.local_key_indices();
      for (size_t ii = 0; ii < coeffs.size(); ++ii) {
        const auto& local_key = coeffs.local_key(ii);
        // find the (sub)entity for this key
        const auto sub_entity = local_key.subEntity();
        const auto codim = local_key.codim();
        const auto& subentity_geometry_type = reference_element.type(sub_entity, codim);
        // and add the respective geometry type
        all_DoF_attached_geometry_types_.insert(subentity_geometry_type);
        geometry_type_to_local_size_[subentity_geometry_type] = local_key_indices[codim][sub_entity].size();
      }
    }
    DUNE_THROW_IF(all_DoF_attached_geometry_types_.size() == 0,
                  Exceptions::mapper_error,
                  "This must not happen, the finite elements report no DoFs attached to (sub)entities!");
    mapper_.update(grid_view_);
  } // ... update_after_adapt(...)

private:
  const GridViewType& grid_view_;
  const LocalFiniteElementFamily& local_finite_elements_;
  const int fe_order_;
  size_t global_size_;
  size_t max_local_size_;
  std::set<GeometryType> all_DoF_attached_geometry_types_;
  std::map<GeometryType, size_t> geometry_type_to_local_size_;
  Implementation mapper_;
}; // class ContinuousMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_CONTINUOUS_HH
