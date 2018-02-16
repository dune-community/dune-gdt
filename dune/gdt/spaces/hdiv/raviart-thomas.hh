// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_HDIV_RAVIART_THOMAS_HH
#define DUNE_GDT_SPACES_HDIV_RAVIART_THOMAS_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/gdt/local/finite-elements/raviart-thomas.hh>
#include <dune/gdt/spaces/basis/raviart-thomas.hh>
#include <dune/gdt/spaces/mapper/continuous.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * The following dimensions/orders/elements are tested to work:
 *
 * - 1d: order 0 works
 * - 2d: order 0 works on simplices, cubes and mixed simplices and cubes
 * - 3d: order 0 work on simplices, cubes
 *
 * The following dimensions/orders/elements are tested to fail:
 *
 * - 3d: mixed simplices and cubes (the mapper cannot handle non-conforming intersections/the switches are not corect)
 */
template <class GV, int p, class R = double>
class RaviartThomasSpace : public SpaceInterface<GV, GV::dimension, 1, R>
{
  static_assert(p == 0, "Higher orders are not tested yet!");
  using ThisType = RaviartThomasSpace<GV, p, R>;
  using BaseType = SpaceInterface<GV, GV::dimension, 1, R>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::GridViewType;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::MapperType;
  using typename BaseType::FiniteElementType;

private:
  using MapperImplementation = ContinuousMapper<GridViewType, FiniteElementType>;
  using GlobalBasisImplementation = RaviartThomasGlobalBasis<GV, R>;

public:
  RaviartThomasSpace(GridViewType grd_vw)
    : grid_view_(grd_vw)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    // create finite elements
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_raviart_thomas_local_finite_element<D, d, R>(geometry_type, p)));
    // check: the mapper does not work for non-conforming intersections
    if (d == 3 && finite_elements_->size() != 1)
      DUNE_THROW(Exceptions::space_error,
                 "when creating a RaviartThomasSpace: non-conforming intersections are not (yet) "
                 "supported, and more than one element type in 3d leads to non-conforming intersections!");
    // compute scaling to ensure basis*integrationElementNormal = 1
    std::map<GeometryType, std::vector<R>> geometry_to_scaling_factors_map;
    for (const auto& geometry_type_and_finite_element_ptr : *finite_elements_) {
      const auto& geometry_type = geometry_type_and_finite_element_ptr.first;
      const auto& finite_element = *geometry_type_and_finite_element_ptr.second;
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto num_intersections = reference_element.size(1);
      geometry_to_scaling_factors_map.insert(std::make_pair(geometry_type, std::vector<R>(num_intersections, R(1.))));
      auto& scaling_factors = geometry_to_scaling_factors_map.at(geometry_type);
      const auto intersection_to_local_DoF_indices_map = finite_element.coefficients().local_key_indices(1);
      for (auto&& intersection_index : XT::Common::value_range(num_intersections)) {
        const auto& normal = reference_element.integrationOuterNormal(intersection_index);
        const auto& intersection_center = reference_element.position(intersection_index, 1);
        const auto shape_functions_evaluations = finite_element.basis().evaluate(intersection_center);
        for (auto&& local_DoF_index : intersection_to_local_DoF_indices_map[intersection_index])
          scaling_factors[intersection_index] /= shape_functions_evaluations[local_DoF_index] * normal;
      }
    }
    // compute switches (as signs of the scaling factor) to ensure continuity of the normal component (therefore we need
    // unique indices for codim 0  entities, which cannot be achieved by the grid layers index set for mixed grids)
    auto entity_indices = std::make_shared<FiniteVolumeMapper<GV>>(grid_view_);
    auto switches = std::make_shared<std::vector<std::vector<R>>>(entity_indices->size());
    for (auto&& entity : elements(grid_view_)) {
      const auto geometry_type = entity.geometry().type();
      const auto& finite_element = *finite_elements_->at(geometry_type);
      const auto& coeffs = finite_element.coefficients();
      const auto entity_index = entity_indices->global_index(entity, 0);
      (*switches)[entity_index] = geometry_to_scaling_factors_map.at(geometry_type);
      auto& local_switches = switches->at(entity_index);
      for (auto&& intersection : intersections(grid_view_, entity)) {
        if (intersection.neighbor() && entity_index < entity_indices->global_index(intersection.outside(), 0)) {
          const auto intersection_index = XT::Common::numeric_cast<unsigned int>(intersection.indexInInside());
          for (size_t ii = 0; ii < coeffs.size(); ++ii) {
            const auto& local_key = coeffs.local_key(ii);
            const auto DoF_subentity_index = local_key.subEntity();
            if (local_key.codim() == 1 && DoF_subentity_index == intersection_index)
              local_switches[DoF_subentity_index] *= -1.;
          }
        }
      }
    }
    // create mapper, basis and communicator
    mapper_ = std::make_shared<MapperImplementation>(grid_view_, finite_elements_);
    basis_ = std::make_shared<GlobalBasisImplementation>(grid_view_, finite_elements_, entity_indices, switches);
    this->create_communicator();
  } // RaviartThomasSpace(...)

  RaviartThomasSpace(const ThisType&) = default;
  RaviartThomasSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const MapperType& mapper() const override final
  {
    return *mapper_;
  }

  const GlobalBasisType& basis() const override final
  {
    return *basis_;
  }

  const FiniteElementType& finite_element(const GeometryType& geometry_type) const override final
  {
    const auto finite_element_search_result = finite_elements_->find(geometry_type);
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen, the grid layer did not report all geometry types!"
                     << "\n   entity.geometry().type() = "
                     << geometry_type);
    return *finite_element_search_result->second;
  }

  SpaceType type() const override final
  {
    return SpaceType::raviart_thomas;
  }

  int min_polorder() const override final
  {
    return p;
  }

  int max_polorder() const override final
  {
    return p;
  }

  bool continuous(const int /*diff_order*/) const override final
  {
    return false;
  }

  bool continuous_normal_components() const override final
  {
    return true;
  }

  bool is_lagrangian() const override final
  {
    return false;
  }


private:
  const GridViewType grid_view_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<MapperImplementation> mapper_;
  std::shared_ptr<GlobalBasisImplementation> basis_;
}; // class RaviartThomasSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_HDIV_RAVIART_THOMAS_HH
