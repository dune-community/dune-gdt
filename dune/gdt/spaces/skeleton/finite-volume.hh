// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_SPACES_SKELETON_FINITE_VOLUME_HH
#define DUNE_GDT_SPACES_SKELETON_FINITE_VOLUME_HH

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basis/finite-volume.hh>
#include <dune/gdt/spaces/mapper/finite-volume-skeleton.hh>
#include <dune/gdt/spaces/interface.hh>


namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class FiniteVolumeSkeletonSpace
{
  static_assert(Dune::AlwaysFalse<GV>::value, "Untested for these dimensions!");
};


template <class GV, class R>
class FiniteVolumeSkeletonSpace<GV, 1, 1, R> : public SpaceInterface<GV, 1, 1, R>
{
  using ThisType = FiniteVolumeSkeletonSpace;
  using BaseType = SpaceInterface<GV, 1, 1, R>;

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
  using MapperImplementation = FiniteVolumeSkeletonMapper<GridViewType, 1, 1>;
  using GlobalBasisImplementation = FiniteVolumeGlobalBasis<GridViewType, 1, R>;

public:
  FiniteVolumeSkeletonSpace(GridViewType grd_vw)
    : BaseType()
    , grid_view_(grd_vw)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, 1>>())
    , mapper_(grid_view_)
    , basis_(grid_view_)
  {
    this->update_after_adapt();
  }

  FiniteVolumeSkeletonSpace(const ThisType& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, 1>>())
    , mapper_(grid_view_)
    , basis_(grid_view_)
  {
    this->update_after_adapt();
  }

  FiniteVolumeSkeletonSpace(ThisType&&) = default;

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
    return mapper_;
  }

  const GlobalBasisType& basis() const override final
  {
    return basis_;
  }

  const LocalFiniteElementFamilyType& finite_elements() const override final
  {
    DUNE_THROW(Exceptions::space_error, "Not implemented yet for skeleton space!");
    return *local_finite_elements_;
  }

  SpaceType type() const override final
  {
    return SpaceType::finite_volume_skeleton;
  }

  int min_polorder() const override final
  {
    return 0;
  }

  int max_polorder() const override final
  {
    return 0;
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

  void restrict_to(
      const ElementType& /*element*/,
      PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& /*persistent_data*/) const override final
  {
    DUNE_THROW(Exceptions::space_error, "Not implemented yet for skeleton space!");
  }

  void update_after_adapt() override final
  {
    // update mapper and basis
    mapper_.update_after_adapt();
    basis_.update_after_adapt();
    this->create_communicator();
  }

  using BaseType::prolong_onto;

  void prolong_onto(const ElementType& /*element*/,
                    const PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& /*persistent_data*/,
                    DynamicVector<R>& /*element_data*/) const override final
  {
    DUNE_THROW(Exceptions::space_error, "Not implemented yet for skeleton space!");
  }

private:
  const GridViewType grid_view_;
  std::unique_ptr<const LocalLagrangeFiniteElementFamily<D, d, R, 1>> local_finite_elements_;
  MapperImplementation mapper_;
  GlobalBasisImplementation basis_;
}; // class FiniteVolumeSkeletonSpace< ..., 1, 1 >


template <size_t r, size_t rC, class R, class GV>
FiniteVolumeSkeletonSpace<GV, r, rC, R> make_finite_volume_skeleton_space(const GV& grid_view)
{
  return FiniteVolumeSkeletonSpace<GV, r, rC, R>(grid_view);
}

template <size_t r, size_t rC, class GV>
FiniteVolumeSkeletonSpace<GV, r, rC, double> make_finite_volume_skeleton_space(const GV& grid_view)
{
  return FiniteVolumeSkeletonSpace<GV, r, rC, double>(grid_view);
}

template <size_t r, class GV>
FiniteVolumeSkeletonSpace<GV, r, 1, double> make_finite_volume_skeleton_space(const GV& grid_view)
{
  return FiniteVolumeSkeletonSpace<GV, r, 1, double>(grid_view);
}

template <class GV>
FiniteVolumeSkeletonSpace<GV, 1, 1, double> make_finite_volume_skeleton_space(const GV& grid_view)
{
  return FiniteVolumeSkeletonSpace<GV, 1, 1, double>(grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_SKELETON_FINITE_VOLUME_HH
