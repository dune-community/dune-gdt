// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
#define DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/type_traits.hh>

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
template <class GV, size_t r = 1, class R = double>
class ContinuousLagrangeSpace : public SpaceInterface<GV, r, 1, R>
{
  using ThisType = ContinuousLagrangeSpace;
  using BaseType = SpaceInterface<GV, r, 1, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalFiniteElementFamilyType;
  using typename BaseType::MapperType;

private:
  using MapperImplementation = ContinuousMapper<GridViewType, LocalFiniteElementFamilyType>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridViewType, r, 1, R>;

public:
  ContinuousLagrangeSpace(GridViewType grd_vw, const int order, const std::string& logging_prefix = "")
    : BaseType(logging_prefix.empty() ? "gdt" : "gdt.spaces.h1.cg",
               logging_prefix.empty() ? "ContinuousLagrangeSpace" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , grid_view_(grd_vw)
    , fe_order_(order)
    , local_finite_elements_(std::make_unique<LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    LOG_(info) << this->logging_id << "(&grd_vw=" << &grd_vw << ", order=" << fe_order_ << ")" << std::endl;
    this->update_after_adapt();
  }

  ContinuousLagrangeSpace(const ThisType& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , fe_order_(other.fe_order_)
    , local_finite_elements_(std::make_unique<LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    this->update_after_adapt();
  }

  ContinuousLagrangeSpace(ThisType&&) = default;

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
    return SpaceType::continuous_lagrange;
  }

  int min_polorder() const override final
  {
    return fe_order_;
  }

  int max_polorder() const override final
  {
    return fe_order_;
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

  void update_after_adapt() override final
  {
    // check: the mapper does not work for non-conforming intersections
    if (d == 3 && grid_view_.indexSet().types(0).size() != 1)
      DUNE_THROW(Exceptions::space_error,
                 "in ContinuousLagrangeSpace: non-conforming intersections are not (yet) "
                 "supported, and more than one element type in 3d leads to non-conforming intersections!");
    // create/update mapper ...
    if (mapper_)
      mapper_->update_after_adapt();
    else
      mapper_ = std::make_unique<MapperImplementation>(grid_view_, *local_finite_elements_, fe_order_);
    // ... and basis
    if (basis_)
      basis_->update_after_adapt();
    else
      basis_ = std::make_unique<GlobalBasisImplementation>(grid_view_, *local_finite_elements_, fe_order_);
    this->create_communicator();
  } // ... update_after_adapt(...)

private:
  const GridViewType grid_view_;
  const int fe_order_;
  int min_polorder_;
  int max_polorder_;
  std::unique_ptr<const LocalLagrangeFiniteElementFamily<D, d, R, r>> local_finite_elements_;
  std::unique_ptr<MapperImplementation> mapper_;
  std::unique_ptr<GlobalBasisImplementation> basis_;
}; // class ContinuousLagrangeSpace


/**
 * \sa ContinuousLagrangeSpace
 */
template <size_t r, class GV, class R = double>
ContinuousLagrangeSpace<GV, r, R> make_continuous_lagrange_space(GV grid_view, const int order)
{
  return ContinuousLagrangeSpace<GV, r, R>(grid_view, order);
}


/**
 * \sa ContinuousLagrangeSpace
 */
template <class GV, class R = double>
ContinuousLagrangeSpace<GV, 1, R> make_continuous_lagrange_space(GV grid_view, const int order)
{
  return ContinuousLagrangeSpace<GV, 1, R>(grid_view, order);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
