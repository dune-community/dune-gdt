// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
#define DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/exceptions.hh>

#include <dune/gdt/local/finite-elements/lagrange.hh>
#include <dune/gdt/spaces/basis/default.hh>
#include <dune/gdt/spaces/mapper/continuous.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * The following dimensions/orders/elements are tested to work:
 *
 * - 1d: orders 1, 2 work
 * - 2d: orders 1, 2 work on simplices, cubes and mixed simplices and cubes
 * - 3d: orders 1, 2 work on simplices, cubes, prisms
 *
 * The following dimensions/orders/elements are tested to fail:
 *
 * - 3d: pyramids (jacobians seem to be incorrect)
 * - 3d: mixed simplices and cubes (intersections are non-conforming)
 *
 * \sa make_local_lagrange_finite_element
 */
template <class GV, class R = double>
class ContinuousLagrangeSpace : public SpaceInterface<GV, 1, 1, R>
{
  using ThisType = ContinuousLagrangeSpace<GV, R>;
  using BaseType = SpaceInterface<GV, 1, 1, R>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::GridViewType;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::MapperType;
  using typename BaseType::FiniteElementType;

private:
  using MapperImplementation = ContinuousMapper<GridViewType, FiniteElementType>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridViewType, 1, 1, R>;

public:
  ContinuousLagrangeSpace(GridViewType grd_vw, const int order)
    : grid_view_(grd_vw)
    , order_(order)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    // create finite elements
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_local_lagrange_finite_element<D, d, R>(geometry_type, order)));
    // check
    if (d == 3 && finite_elements_->size() != 1)
      DUNE_THROW(Exceptions::space_error,
                 "when creating a ContinuousLagrangeSpace: non-conforming intersections are not (yet) "
                 "supported, and more than one element type in 3d leads to non-conforming intersections!");
    // create mapper, basis and communicator
    mapper_ = std::make_shared<MapperImplementation>(grid_view_, finite_elements_);
    basis_ = std::make_shared<GlobalBasisImplementation>(grid_view_, finite_elements_);
    this->create_communicator();
  } // ContinuousLagrangeSpace(...)

  ContinuousLagrangeSpace(const ThisType&) = default;
  ContinuousLagrangeSpace(ThisType&&) = default;

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
    return SpaceType::continuous_lagrange;
  }

  int min_polorder() const override final
  {
    return order_;
  }

  int max_polorder() const override final
  {
    return order_;
  }

  bool continuous(const int diff_order) const override final
  {
    return diff_order == 0;
  }

  bool continuous_normal_components() const override final
  {
    return true;
  }

  bool is_lagrangian() const override final
  {
    return true;
  }

private:
  const GridViewType grid_view_;
  const int order_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<MapperImplementation> mapper_;
  std::shared_ptr<GlobalBasisImplementation> basis_;
}; // class ContinuousLagrangeSpace


/**
 * \sa ContinuousLagrangeSpace
 */
template <class GV, class R = double>
ContinuousLagrangeSpace<GridView<GV>, R> make_continuous_lagrange_space(GridView<GV> grid_view, const int order)
{
  return ContinuousLagrangeSpace<GridView<GV>, R>(grid_view, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
