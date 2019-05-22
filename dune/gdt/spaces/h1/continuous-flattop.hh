// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_SPACES_H1_CONTINUOUS_FLATTOP_HH
#define DUNE_GDT_SPACES_H1_CONTINUOUS_FLATTOP_HH

#include <memory>
//#include <vector>

//#include <dune/common/typetraits.hh>

//#include <dune/geometry/type.hh>

#include <dune/grid/common/gridview.hh>

//#include <dune/xt/common/exceptions.hh>
//#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/finite-elements/flattop.hh>
#include <dune/gdt/spaces/basis/default.hh>
#include <dune/gdt/spaces/mapper/continuous.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * \sa make_local_lagrange_finite_element
 */
template <class GV, size_t r = 1, class R = double>
class ContinuousFlatTopSpace : public SpaceInterface<GV, r, 1, R>
{
  using ThisType = ContinuousFlatTopSpace;
  using BaseType = SpaceInterface<GV, r, 1, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::GridViewType;
  using typename BaseType::LocalFiniteElementFamilyType;
  using typename BaseType::MapperType;

private:
  using MapperImplementation = ContinuousMapper<GridViewType, LocalFiniteElementFamilyType, r>;
  using GlobalBasisImplementation = DefaultGlobalBasis<GridViewType, r, 1, R>;

public:
  ContinuousFlatTopSpace(GridViewType grd_vw, const int fe_order, const D& overlap = 0.5)
    : grid_view_(grd_vw)
    , fe_order_(fe_order)
    , local_finite_elements_(std::make_unique<LocalFlatTopFiniteElementFamily<D, d, R, r>>(overlap))
    , mapper_(nullptr)
    , basis_(nullptr)
  {
    this->update_after_adapt();
  }

  ContinuousFlatTopSpace(const ThisType&) = default;
  ContinuousFlatTopSpace(ThisType&&) = default;

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
    return SpaceType::continuous_lagrange;
  }

  int min_polorder() const override final
  {
    return fe_order_;
  }

  int max_polorder() const override final
  {
    return fe_order_ + 1;
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
                 "in ContinuousFlatTopSpace: non-conforming intersections are not (yet) "
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
  std::unique_ptr<const LocalFlatTopFiniteElementFamily<D, d, R, r>> local_finite_elements_;
  std::unique_ptr<MapperImplementation> mapper_;
  std::unique_ptr<GlobalBasisImplementation> basis_;
}; // class ContinuousFlatTopSpace


/**
 * \sa ContinuousFlatTopSpace
 */
template <size_t r, class GV, class R = double>
ContinuousFlatTopSpace<GV, r, R>
make_continuous_flattop_space(GV grid_view, const int order, const double& overlap = 0.5)
{
  return ContinuousFlatTopSpace<GV, r, R>(grid_view, order, overlap);
}


/**
 * \sa ContinuousFlatTopSpace
 */
template <class GV, class R = double>
ContinuousFlatTopSpace<GV, 1, R>
make_continuous_flattop_space(GV grid_view, const int order, const double& overlap = 0.5)
{
  return ContinuousFlatTopSpace<GV, 1, R>(grid_view, order, overlap);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_H1_CONTINUOUS_FLATTOP_HH
