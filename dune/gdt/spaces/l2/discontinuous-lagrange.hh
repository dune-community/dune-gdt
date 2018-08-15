// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH
#define DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>

#include <dune/gdt/local/finite-elements/lagrange.hh>
#include <dune/gdt/spaces/basis/default.hh>
#include <dune/gdt/spaces/mapper/discontinuous.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * The following dimensions/orders/elements are tested to work:
 *
 * - 1d: orders 0, ..., 18
 * - 2d: orders 0, ..., 10 work on simplices, cubes and mixed simplices and cubes
 * - 2d: orders 11, ..., 15 also work on simplices
 * - 3d: orders 0, ..., 7 work on simplices, cubes, prisms and mixed simplices and cubes
 * - 3d: orders 8, ..., 14 also work on simplices
 * - 3d: orders 8, 9 also work on prisms
 *
 * The following dimensions/orders/elements are tested to fail (basis matrix fails to invert):
 *
 * - 1d: orders > 18
 * - 2d: orders > 15 on simplices
 * - 2d: orders > 10 on cubes
 * - 3d: orders > 14 on simplices
 * - 3d: orders > 7 on cubes
 * - 3d: orders > 9 on prisms
 *
 * \sa make_local_lagrange_finite_element
 * \sa make_discontinuous_lagrange_space
 */
template <class GV, size_t r = 1, class R = double>
class DiscontinuousLagrangeSpace : public SpaceInterface<GV, r, 1, R>
{
  using ThisType = DiscontinuousLagrangeSpace<GV, r, R>;
  using BaseType = SpaceInterface<GV, r, 1, R>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::GridViewType;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::MapperType;
  using typename BaseType::FiniteElementType;

private:
  using MapperImplementation = DiscontinuousMapper<GridViewType, FiniteElementType>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridViewType, r, 1, R>;

public:
  DiscontinuousLagrangeSpace(GridViewType grd_vw, const int order)
    : grid_view_(grd_vw)
    , order_(order)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    // create finite elements
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_local_lagrange_finite_element<D, d, R, r>(geometry_type, order_)));
    // create mapper, basis and communicator
    mapper_ = std::make_shared<MapperImplementation>(grid_view_, finite_elements_);
    basis_ = std::make_shared<GlobalBasisImplementation>(grid_view_, finite_elements_);
    this->create_communicator();
  }

  DiscontinuousLagrangeSpace(const ThisType&) = default;
  DiscontinuousLagrangeSpace(ThisType&&) = default;

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
    return SpaceType::discontinuous_lagrange;
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
    return false;
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
}; // class DiscontinuousLagrangeSpace


/**
 * \sa DiscontinuousLagrangeSpace
 */
template <class GV, class R = double>
DiscontinuousLagrangeSpace<GridView<GV>, 1, R> make_discontinuous_lagrange_space(GridView<GV> grid_view,
                                                                                 const int order)
{
  return DiscontinuousLagrangeSpace<GridView<GV>, 1, R>(grid_view, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH
