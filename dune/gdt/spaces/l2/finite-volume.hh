// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

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
  using ThisType = FiniteVolumeSpace<GV, r, 1, R>;
  using BaseType = SpaceInterface<GV, r, 1, R>;

public:
  using typename BaseType::D;
  using BaseType::d;
  using typename BaseType::G;
  using typename BaseType::ElementType;
  using typename BaseType::GridViewType;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::MapperType;
  using typename BaseType::FiniteElementType;

private:
  using MapperImplementation = FiniteVolumeMapper<GridViewType, r, 1>;
  using GlobalBasisImplementation = FiniteVolumeGlobalBasis<GridViewType, r, R>;

public:
  FiniteVolumeSpace(GridViewType grd_vw)
    : grid_view_(grd_vw)
    , finite_elements_(new std::map<GeometryType, std::shared_ptr<FiniteElementType>>())
    , mapper_(new MapperImplementation(grid_view_))
    , basis_(new GlobalBasisImplementation(grid_view_))
  {
    // create finite elements
    for (auto&& geometry_type : grid_view_.indexSet().types(0))
      finite_elements_->insert(
          std::make_pair(geometry_type, make_local_lagrange_finite_element<D, d, R, r>(geometry_type, 0)));
    // create communicator
    this->create_communicator();
  }

  FiniteVolumeSpace(const ThisType&) = default;
  FiniteVolumeSpace(ThisType&&) = default;

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
   * On the domain covered by the coarser element, this computes an L^2 projection of the function defined on the finer
   * elements (which corresponds to weighted averaging in the FV case).
   */
  void restrict_to(const ElementType& element,
                   PersistentContainer<G, DynamicVector<R>>& persistent_data) const override final
  {
    auto& element_restriction_data = persistent_data[element];
    if (element_restriction_data.size() == 0) {
      DUNE_THROW_IF(element.isLeaf(), Exceptions::space_error, "");
      for (auto&& child_element : descendantElements(element, element.level() + 1)) {
        // ensure we have data on all descendant elements of the next level
        this->restrict_to(child_element, persistent_data);
        // compute restriction
        auto child_restriction_data = persistent_data[child_element];
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

protected:
  /**
   * \note In general, we would have to check for newly created GeometryTypes and to recreate the local FEs accordingly.
   *       This is postponed until we have the LocalFiniteElementFamily.
   */
  void update_after_adapt() override final
  {
    basis_->update_after_adapt();
    mapper_->update_after_adapt();
  }

public:
  using BaseType::prolong_onto;

  void prolong_onto(const ElementType& element,
                    const PersistentContainer<G, DynamicVector<R>>& persistent_data,
                    DynamicVector<R>& element_data) const override final
  {
    if (element.isNew()) {
      // ... by prolongation from the father element ...
      const auto& father_data = persistent_data[element.father()];
      DUNE_THROW_IF(father_data.size() == 0, Exceptions::space_error, "");
      if (element_data.size() != father_data.size())
        element_data.resize(father_data.size());
      element_data = father_data;
    } else {
      // ... or by copying the data from this unchanged element
      const auto& original_element_data = persistent_data[element];
      DUNE_THROW_IF(original_element_data.size() == 0, Exceptions::space_error, "");
      if (element_data.size() != original_element_data.size())
        element_data.resize(original_element_data.size());
      element_data = original_element_data;
    }
  } // ... prolong_onto(...)

private:
  const GridViewType grid_view_;
  std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  std::shared_ptr<MapperImplementation> mapper_;
  const std::shared_ptr<GlobalBasisImplementation> basis_;
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
