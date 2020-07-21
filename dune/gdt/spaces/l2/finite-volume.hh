// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014, 2016 - 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016, 2018)

#ifndef DUNE_GDT_SPACES_L2_FINITE_VOLUME_HH
#define DUNE_GDT_SPACES_L2_FINITE_VOLUME_HH

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basis/finite-volume.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>
#include <dune/gdt/spaces/interface.hh>


namespace Dune {
namespace GDT {


// forward, to allow for specialization
template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class FiniteVolumeSpace
{
  static_assert(Dune::AlwaysFalse<GV>::value, "Untested for these dimensions!");
};


template <class GV, size_t r, class R>
class FiniteVolumeSpace<GV, r, 1, R> : public SpaceInterface<GV, r, 1, R>
{
  using ThisType = FiniteVolumeSpace;
  using BaseType = SpaceInterface<GV, r, 1, R>;

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
  using MapperImplementation = FiniteVolumeMapper<GridViewType, r, 1>;
  using GlobalBasisImplementation = FiniteVolumeGlobalBasis<GridViewType, r, R>;

public:
  FiniteVolumeSpace(GridViewType grd_vw)
    : BaseType()
    , grid_view_(grd_vw)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(grid_view_)
    , basis_(grid_view_)
  {
    this->update_after_adapt();
  }

  FiniteVolumeSpace(const ThisType& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , local_finite_elements_(std::make_unique<const LocalLagrangeFiniteElementFamily<D, d, R, r>>())
    , mapper_(grid_view_)
    , basis_(grid_view_)
  {
    this->update_after_adapt();
  }

  FiniteVolumeSpace(ThisType&&) = default;

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
    return *local_finite_elements_;
  }

  SpaceType type() const override final
  {
    return SpaceType::finite_volume;
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

  /**
   * More efficient restriction than in SpaceInterface for FV.
   */
  void restrict_to(
      const ElementType& element,
      PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& persistent_data) const override final
  {
    auto& element_restriction_data = persistent_data[element].second;
    if (element_restriction_data.size() == 0) {
      DUNE_THROW_IF(element.isLeaf(), Exceptions::space_error, "");
      for (auto&& child_element : descendantElements(element, element.level() + 1)) {
        // ensure we have data on all descendant elements of the next level
        this->restrict_to(child_element, persistent_data);
        // compute restriction
        auto child_restriction_data = persistent_data[child_element].second;
        child_restriction_data *= child_element.geometry().volume();
        if (element_restriction_data.size() == 0)
          element_restriction_data = child_restriction_data; // initialize with child data
        else
          element_restriction_data += child_restriction_data;
      }
      // now we have assembled h*value
      element_restriction_data /= element.geometry().volume();
    }
  } // ... restrict_to(...)

  void update_after_adapt() override final
  {
    // update mapper and basis
    mapper_.update_after_adapt();
    basis_.update_after_adapt();
    this->create_communicator();
  }

  using BaseType::prolong_onto;

  /**
   * More efficient prolongation than in SpaceInterface for FV.
   */
  void prolong_onto(const ElementType& element,
                    const PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& persistent_data,
                    DynamicVector<R>& element_data) const override final
  {
    if (element.isNew()) {
      // ... by prolongation from the father element ...
      const auto& father_data = persistent_data[element.father()].second;
      DUNE_THROW_IF(father_data.size() == 0, Exceptions::space_error, "");
      if (element_data.size() != father_data.size())
        element_data.resize(father_data.size());
      element_data = father_data;
    } else {
      // ... or by copying the data from this unchanged element
      const auto& original_element_data = persistent_data[element].second;
      DUNE_THROW_IF(original_element_data.size() == 0, Exceptions::space_error, "");
      if (element_data.size() != original_element_data.size())
        element_data.resize(original_element_data.size());
      element_data = original_element_data;
    }
  } // ... prolong_onto(...)

private:
  const GridViewType grid_view_;
  std::unique_ptr<const LocalLagrangeFiniteElementFamily<D, d, R, r>> local_finite_elements_;
  MapperImplementation mapper_;
  GlobalBasisImplementation basis_;
}; // class FiniteVolumeSpace< ..., r, 1 >


template <size_t r, size_t rC, class R, class GV>
FiniteVolumeSpace<GV, r, rC, R> make_finite_volume_space(const GV& grid_view)
{
  return FiniteVolumeSpace<GV, r, rC, R>(grid_view);
}

template <size_t r, size_t rC, class GV>
FiniteVolumeSpace<GV, r, rC, double> make_finite_volume_space(const GV& grid_view)
{
  return FiniteVolumeSpace<GV, r, rC, double>(grid_view);
}

template <size_t r, class GV>
FiniteVolumeSpace<GV, r, 1, double> make_finite_volume_space(const GV& grid_view)
{
  return FiniteVolumeSpace<GV, r, 1, double>(grid_view);
}

template <class GV>
FiniteVolumeSpace<GV, 1, 1, double> make_finite_volume_space(const GV& grid_view)
{
  return FiniteVolumeSpace<GV, 1, 1, double>(grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_L2_FINITE_VOLUME_HH
