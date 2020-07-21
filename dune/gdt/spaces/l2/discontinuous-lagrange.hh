// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

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
  using ThisType = DiscontinuousLagrangeSpace;
  using BaseType = SpaceInterface<GV, r, 1, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalFiniteElementFamilyType;
  using typename BaseType::MapperType;

private:
  using MapperImplementation = DiscontinuousMapper<GridViewType, LocalFiniteElementFamilyType>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridViewType, r, 1, R>;

public:
  DiscontinuousLagrangeSpace(GridViewType grd_vw, const int order = 1)
    : BaseType()
    , grid_view_(grd_vw)
    , order_(order)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    this->update_after_adapt();
  }

  DiscontinuousLagrangeSpace(const ThisType& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , order_(other.order_)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    this->update_after_adapt();
  }

  DiscontinuousLagrangeSpace(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;

  ThisType& operator=(ThisType&&) = delete;

  BaseType* copy() const override final
  {
    return new ThisType(*this);
  }

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

  void update_after_adapt() override final
  {
    // create/update mapper ...
    if (mapper_)
      mapper_->update_after_adapt();
    else
      mapper_ = std::make_unique<MapperImplementation>(grid_view_, *local_finite_elements_, order_);
    // ... and basis
    if (basis_)
      basis_->update_after_adapt();
    else
      basis_ = std::make_unique<GlobalBasisImplementation>(grid_view_, *local_finite_elements_, order_);
    this->create_communicator();
  } // ... update_after_adapt(...)

private:
  const GridViewType grid_view_;
  const int order_;
  std::unique_ptr<const LocalLagrangeFiniteElementFamily<D, d, R, r>> local_finite_elements_;
  std::unique_ptr<MapperImplementation> mapper_;
  std::unique_ptr<GlobalBasisImplementation> basis_;
}; // class DiscontinuousLagrangeSpace


/**
 * \sa DiscontinuousLagrangeSpace
 */
template <size_t r, class GV, class R = double>
DiscontinuousLagrangeSpace<GV, r, R> make_discontinuous_lagrange_space(GV grid_view, const int order)
{
  return DiscontinuousLagrangeSpace<GV, r, R>(grid_view, order);
}


/**
 * \sa DiscontinuousLagrangeSpace
 */
template <class GV, class R = double>
DiscontinuousLagrangeSpace<GV, 1, R> make_discontinuous_lagrange_space(GV grid_view, const int order)
{
  return DiscontinuousLagrangeSpace<GV, 1, R>(grid_view, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH
