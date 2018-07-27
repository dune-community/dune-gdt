// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_HH

#include <mutex>

#include <dune/common/unused.hh>

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/crtp.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/la/container/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace internal {


//// forward, needed for friendlyness
// template <class TestSpaceType, class AnsatzSpaceType, class GridLayerType, class ConstraintsType>
// class ConstraintsWrapper;


} // namespace internal


/**
 * \brief CRTP interface for all implementations of constraints.
 *
 *        We need this interface for template matching in the SystemAssembler.
 */
// template <class Traits>
// class ConstraintsInterface : public XT::CRTPInterface<ConstraintsInterface<Traits>, Traits>
//{
// public:
//  typedef typename Traits::derived_type derived_type;
//}; // class ConstraintsInterface


//// forward
// template <class IntersectionType>
// class DirichletConstraints;


// namespace internal {


// template <class IntersectionType>
// class DirichletConstraintsTraits
//{
// public:
//  typedef DirichletConstraints<IntersectionType> derived_type;
//};


//} // namespace internal


template <class IntersectionType, class SpaceType>
class DirichletConstraints : public Dune::XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using ThisType = DirichletConstraints<IntersectionType, SpaceType>;
  using BaseType = XT::Grid::ElementFunctor<typename SpaceType::GridViewType>;
  //  using Propagator = XT::Common::ThreadResultPropagator<ThisType, size_t>;
  //  friend Propagator;


public:
  using BoundaryInfoType = XT::Grid::BoundaryInfo<IntersectionType>;
  using ElementType = typename ThisType::ElementType;
  using GridView = typename SpaceType::GridViewType;
  static const constexpr size_t r = SpaceType::r;
  static const constexpr size_t rC = SpaceType::rC;
  using R = typename SpaceType::R;

  DirichletConstraints(const BoundaryInfoType& bnd_info, const SpaceType& space, const bool set = true)
    : boundary_info_(bnd_info)
    , space_(space)
    , set_(set)
  {
  }

  void apply_local(const ElementType& element) override final
  {
    std::set<size_t> local_DoFs;
    static const XT::Grid::DirichletBoundary dirichlet{};
    const auto lps = space_.finite_element(element.type()).lagrange_points();
    const auto intersection_it_end = space_.grid_view().iend(element);
    for (auto intersection_it = space_.grid_view().ibegin(element); intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundary_info_.type(intersection) == dirichlet || (!intersection.neighbor() && !intersection.boundary()))
        for (size_t ii = 0; ii < lps.size(); ++ii) {
          const auto& local_lagrange_point = lps[ii];
          if (XT::Grid::contains(intersection, element.geometry().global(local_lagrange_point)))
            local_DoFs.insert(ii);
        }
    }

    if (local_DoFs.size() > 0) {
      for (const auto& local_DoF : local_DoFs) {
        dirichlet_DoFs_.insert(space_.mapper().global_index(element, local_DoF));
      }
    }
  }

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  //  size_t size() const
  //  {
  //    return size_;
  //  }

  //  inline void insert(const size_t DoF)
  //  {
  //    DUNE_UNUSED std::lock_guard<std::mutex> mutex_guard(mutex_);
  //    assert(DoF < size_);
  //    dirichlet_DoFs_.insert(DoF);
  //  }

  const std::set<size_t>& dirichlet_DoFs() const
  {
    return dirichlet_DoFs_;
  }

  template <class M>
  void apply(XT::LA::MatrixInterface<M>& matrix) const
  {
    if (set_) {
      for (const auto& DoF : dirichlet_DoFs_)
        matrix.unit_row(DoF);
    } else {
      for (const auto& DoF : dirichlet_DoFs_)
        matrix.clear_row(DoF);
    }
  } // ... apply(...)

  template <class V>
  void apply(XT::LA::VectorInterface<V>& vector) const
  {
    for (const auto& DoF : dirichlet_DoFs_)
      vector[DoF] = 0.0;
  }

  template <class M, class V>
  void apply(XT::LA::MatrixInterface<M>& matrix, XT::LA::VectorInterface<V>& vector) const
  {
    if (set_) {
      for (const auto& DoF : dirichlet_DoFs_) {
        matrix.unit_row(DoF);
        vector[DoF] = 0.0;
      }
    } else {
      for (const auto& DoF : dirichlet_DoFs_) {
        matrix.clear_row(DoF);
        vector[DoF] = 0.0;
      }
    }
  } // ... apply(...)

  //  void finalize() override
  //  {
  //    Propagator::finalize_imp();
  //  }

  BaseType* copy() override
  {
    return new ThisType(*this);
    //    return Propagator::copy_imp();
  }

private:
  const BoundaryInfoType& boundary_info_;
  const SpaceInterface<GridView, r, rC, R>& space_;
  //  const size_t size_;
  const bool set_;
  std::set<size_t> dirichlet_DoFs_;
  //  std::mutex mutex_;
}; // class DirichletConstraints


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
