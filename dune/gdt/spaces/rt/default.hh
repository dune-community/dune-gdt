// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_RT_DEFAULT_HH
#define DUNE_GDT_SPACES_RT_DEFAULT_HH

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


// forwards, required for the traits
template <class GL, int p, class R = double>
class RaviartThomasSpace;


namespace internal {


template <class GL, int p, class R>
class RaviartThomasSpaceTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(p == 0, "Not implemented yet!");
  using G = XT::Grid::extract_grid_t<GL>;

public:
  static const constexpr int polOrder = p;
  static const constexpr size_t dimDomain = GL::dimension;
  static const constexpr size_t dimRange = dimDomain;
  static const constexpr size_t dimRangeCols = 1;

  using derived_type = RaviartThomasSpace<GL, p, R>;
  static const constexpr bool continuous = false;
  using GridLayerType = GL;
  using BackendType = double;
  using RangeFieldType = R;
  static const constexpr XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const constexpr Backends backend_type{Backends::gdt};
  typedef DofCommunicationChooser<GridLayerType, false> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;

private:
  friend class RaviartThomasSpace<GL, p, R>;
}; // class RaviartThomasSpaceTraits


} // namespace internal


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
template <class GL, int p, class R>
class RaviartThomasSpace
    : public SpaceInterface<internal::RaviartThomasSpaceTraits<GL, p, R>, GL::dimension, GL::dimension>
{
public:
  using Traits = internal::RaviartThomasSpaceTraits<GL, p, R>;

private:
  using BaseType = SpaceInterface<internal::RaviartThomasSpaceTraits<GL, p, R>, GL::dimension, GL::dimension>;
  using ThisType = RaviartThomasSpace<GL, p, R>;
  using D = typename GL::ctype;
  static const constexpr size_t d = BaseType::dimDomain;

public:
  using typename BaseType::FiniteElementType;
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::MapperType;
  using typename BaseType::GlobalBasisType;
  typedef typename Traits::DofCommunicatorType DofCommunicatorType;

private:
  typedef typename Traits::DofCommunicationChooserType DofCommunicationChooserType;
  using MapperImplementation = ContinuousMapper<GL, FiniteElementType>;
  using GlobalBasisImplementation = RaviartThomasGlobalBasis<GL, R>;

public:
  RaviartThomasSpace(GridLayerType grd_lr)
    : grid_layer_(grd_lr)
    , communicator_(DofCommunicationChooserType::create(grid_layer_))
    , backend_(0)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , geometry_to_local_DoF_indices_map_(new std::map<GeometryType, std::vector<size_t>>())
    , mapper_(nullptr)
    , basis_(nullptr)
    , communicator_prepared_(false)
  {
    // create finite elements
    for (auto&& geometry_type : grid_layer_.indexSet().types(0)) {
      auto finite_element = make_raviart_thomas_local_finite_element<D, d, R>(geometry_type, p);
      // this is only valid for p0 elements
      geometry_to_local_DoF_indices_map_->insert(
          std::make_pair(geometry_type, std::vector<size_t>(finite_element->size())));
      finite_elements_->insert(std::make_pair(geometry_type, std::move(finite_element)));
    }
    // check: the mapper does not work for non-conforming intersections
    if (d == 3 && finite_elements_->size() != 1)
      DUNE_THROW(space_error,
                 "when creating a RaviartThomasSpace: non-conforming intersections are not (yet) "
                 "supported, and more than one element type in 3d leads to non-conforming intersections!");
    // compute local-key-to-intersection relationship
    for (const auto& geometry_type_and_finite_element_ptr : *finite_elements_) {
      const auto& geometry_type = geometry_type_and_finite_element_ptr.first;
      const auto& finite_element = *geometry_type_and_finite_element_ptr.second;
      const auto& coeffs = finite_element.coefficients();
      auto& local_key_to_intersection_map = geometry_to_local_DoF_indices_map_->at(geometry_type);
      for (size_t ii = 0; ii < coeffs.size(); ++ii) {
        const auto& local_key = coeffs.local_key(ii);
        if (local_key.index() != 0)
          DUNE_THROW(XT::Common::Exceptions::internal_error, "This must not happen for p0!");
        if (local_key.codim() != 1)
          DUNE_THROW(XT::Common::Exceptions::internal_error, "This must not happen for p0!");
        local_key_to_intersection_map[ii] = local_key.subEntity();
      }
    }
    // compute scaling to ensure basis*integrationElementNormal = 1
    std::map<GeometryType, std::vector<R>> geometry_to_scaling_factors_map;
    for (const auto& geometry_type_and_finite_element_ptr : *finite_elements_) {
      const auto& geometry_type = geometry_type_and_finite_element_ptr.first;
      const auto& finite_element = *geometry_type_and_finite_element_ptr.second;
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto num_intersections = reference_element.size(1);
      geometry_to_scaling_factors_map.insert(std::make_pair(geometry_type, std::vector<R>(num_intersections, R(1.))));
      auto& scaling_factors = geometry_to_scaling_factors_map.at(geometry_type);
      const auto& DoF_indices = geometry_to_local_DoF_indices_map_->at(geometry_type);
      for (auto&& intersection_index : XT::Common::value_range(num_intersections)) {
        const auto& normal = reference_element.integrationOuterNormal(intersection_index);
        const auto& intersection_center = reference_element.position(intersection_index, 1);
        const auto DoF_index = DoF_indices[intersection_index];
        const auto shape_functions_evaluations = finite_element.basis().evaluate(intersection_center);
        scaling_factors[intersection_index] /= shape_functions_evaluations[DoF_index] * normal;
      }
    }
    // compute switches (as signs of the scaling factor) to ensure continuity of the normal component (therefore we need
    // unique indices for codim 0  entities, which cannot be achieved by the grid layers index set for mixed grids)
    auto entity_indices = std::make_shared<FiniteVolumeMapper<GL>>(grid_layer_);
    auto switches = std::make_shared<std::vector<std::vector<R>>>(entity_indices->size());
    for (auto&& entity : elements(grid_layer_)) {
      const auto geometry_type = entity.geometry().type();
      const auto& finite_element = *finite_elements_->at(geometry_type);
      const auto& coeffs = finite_element.coefficients();
      const auto entity_index = entity_indices->global_index(entity, 0);
      (*switches)[entity_index] = geometry_to_scaling_factors_map.at(geometry_type);
      auto& local_switches = switches->at(entity_index);
      for (auto&& intersection : intersections(grid_layer_, entity)) {
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
    // create mapper and basis
    mapper_ = std::make_shared<MapperImplementation>(grid_layer_, finite_elements_);
    basis_ = std::make_shared<GlobalBasisImplementation>(finite_elements_, entity_indices, switches);
  } // RaviartThomasSpace(...)

  RaviartThomasSpace(const ThisType&) = default;
  RaviartThomasSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  const GlobalBasisType& basis() const
  {
    return *basis_;
  }

  DofCommunicatorType& dof_communicator() const
  {
    if (!communicator_prepared_) {
      communicator_prepared_ = DofCommunicationChooserType::prepare(*this, *communicator_);
    }
    return *communicator_;
  }

  /// \note this makes sense only for for p0, once we export the local finite element, this should become obsolete!
  std::vector<size_t> local_DoF_indices(const EntityType& entity) const
  {
    const auto search_result = geometry_to_local_DoF_indices_map_->find(entity.geometry().type());
    if (search_result == geometry_to_local_DoF_indices_map_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen after the checks in the ctor, the grid layer did not report all geometry types!"
                 "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    return search_result->second;
  } // ... local_DoF_indices(...)

private:
  const GridLayerType grid_layer_;
  mutable std::shared_ptr<DofCommunicatorType> communicator_;
  const double backend_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<std::map<GeometryType, std::vector<size_t>>> geometry_to_local_DoF_indices_map_;
  std::shared_ptr<MapperImplementation> mapper_;
  std::shared_ptr<GlobalBasisImplementation> basis_;
  mutable bool communicator_prepared_;
}; // class RaviartThomasSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_DEFAULT_HH
