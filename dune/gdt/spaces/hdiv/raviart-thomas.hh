// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

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
template <class GV, class R = double>
class RaviartThomasSpace : public SpaceInterface<GV, GV::dimension, 1, R>
{
  using ThisType = RaviartThomasSpace;
  using BaseType = SpaceInterface<GV, GV::dimension, 1, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::ElementType;
  using typename BaseType::G;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalFiniteElementFamilyType;
  using typename BaseType::MapperType;

private:
  using MapperImplementation = ContinuousMapper<GV, LocalFiniteElementFamilyType>;
  using GlobalBasisImplementation = RaviartThomasGlobalBasis<GV, R>;

public:
  RaviartThomasSpace(GridViewType grd_vw, const int order)
    : grid_view_(grd_vw)
    , order_(order)
    , local_finite_elements_(std::make_unique<const LocalRaviartThomasFiniteElementFamily<D, d, R>>())
    , element_indices_(grid_view_)
    , fe_data_()
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    DUNE_THROW_IF(order_ != 0, Exceptions::space_error, "Higher orders are not testet yet!");
    this->update_after_adapt();
  }

  RaviartThomasSpace(const ThisType& other)
    : grid_view_(other.grid_view_)
    , order_(other.order_)
    , local_finite_elements_(std::make_unique<const LocalRaviartThomasFiniteElementFamily<D, d, R>>())
    , element_indices_(grid_view_)
    , fe_data_()
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    this->update_after_adapt();
  }

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

  RaviartThomasSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const GridViewType& grid_view() const override final
  {
    return grid_view_;
  }

  const MapperType& mapper() const override final
  {
    assert(mapper_ && "This must not happen!");
    return *mapper_;
  }

  const GlobalBasisType& basis() const override final
  {
    assert(basis_ && "This must not happen!");
    return *basis_;
  }

  const LocalFiniteElementFamilyType& finite_elements() const override final
  {
    return *local_finite_elements_;
  }

  SpaceType type() const override final
  {
    return SpaceType::raviart_thomas;
  }

  int min_polorder() const override final
  {
    return order_;
  }

  int max_polorder() const override final
  {
    return order_;
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

protected:
  void update_after_adapt() override final
  {
    element_indices_.update_after_adapt();
    // check: the mapper does not work for non-conforming intersections
    if (d == 3 && grid_view_.indexSet().types(0).size() != 1)
      DUNE_THROW(Exceptions::space_error,
                 "in RaviartThomasSpace: non-conforming intersections are not (yet) "
                 "supported, and more than one element type in 3d leads to non-conforming intersections!");
    // compute scaling to ensure basis*integrationElementNormal = 1
    std::map<GeometryType, DynamicVector<R>> geometry_to_scaling_factors_map;
    for (const auto& geometry_type : grid_view_.indexSet().types(0)) {
      const auto& finite_element = local_finite_elements_->get(geometry_type, order_);
      const auto& reference_element = ReferenceElements<D, d>::general(geometry_type);
      const auto num_intersections = reference_element.size(1);
      geometry_to_scaling_factors_map.insert(std::make_pair(geometry_type, DynamicVector<R>(num_intersections, R(1.))));
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
    fe_data_.resize(element_indices_.size());
    for (auto&& entity : elements(grid_view_)) {
      const auto geometry_type = entity.type();
      const auto& finite_element = local_finite_elements_->get(geometry_type, order_);
      const auto& coeffs = finite_element.coefficients();
      const auto element_index = element_indices_.global_index(entity, 0);
      fe_data_[element_index] = geometry_to_scaling_factors_map.at(geometry_type);
      auto& local_switches = fe_data_[element_index];
      for (auto&& intersection : intersections(grid_view_, entity)) {
        if (intersection.neighbor() && element_index < element_indices_.global_index(intersection.outside(), 0)) {
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
    // create/update mapper ...
    if (mapper_)
      mapper_->update_after_adapt();
    else
      mapper_ = std::make_unique<MapperImplementation>(grid_view_, *local_finite_elements_, order_);
    // ... and basis
    if (basis_)
      basis_->update_after_adapt();
    else
      basis_ = std::make_unique<GlobalBasisImplementation>(
          grid_view_, order_, *local_finite_elements_, element_indices_, fe_data_);
    this->create_communicator();
  } // ... update_after_adapt(...)

private:
  const GridViewType grid_view_;
  const int order_;
  std::unique_ptr<const LocalRaviartThomasFiniteElementFamily<D, d, R>> local_finite_elements_;
  FiniteVolumeMapper<GV> element_indices_;
  std::vector<DynamicVector<R>> fe_data_;
  std::unique_ptr<MapperImplementation> mapper_;
  std::unique_ptr<GlobalBasisImplementation> basis_;
}; // class RaviartThomasSpace


/**
 * \sa RaviartThomasSpace
 */
template <class GV, class R = double>
RaviartThomasSpace<GV, R> make_raviart_thomas_space(GV grid_view, const int order)
{
  return RaviartThomasSpace<GV, R>(grid_view, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_HDIV_RAVIART_THOMAS_HH
