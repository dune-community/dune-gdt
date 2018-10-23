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


template <typename T>
struct setUnion
{
  std::set<T> operator()(const std::set<T>& a, const std::set<T>& b)
  {
    std::set<T> result = a;
    result.insert(b.begin(), b.end());
    return result;
  }
};


} // namespace internal


template <class IntersectionType, class SpaceType>
class DirichletConstraints
    : public Dune::XT::Grid::ElementFunctor<typename SpaceType::GridViewType>,
      public XT::Common::ThreadResultPropagator<DirichletConstraints<IntersectionType, SpaceType>,
                                                std::set<size_t>,
                                                internal::setUnion<size_t>>
{
  using ThisType = DirichletConstraints<IntersectionType, SpaceType>;
  using BaseType = XT::Grid::ElementFunctor<typename SpaceType::GridViewType>;
  using Propagator = XT::Common::ThreadResultPropagator<ThisType, std::set<size_t>, internal::setUnion<size_t>>;
  friend Propagator;

public:
  using BoundaryInfoType = XT::Grid::BoundaryInfo<IntersectionType>;
  using ElementType = typename ThisType::ElementType;
  using GridView = typename SpaceType::GridViewType;
  static const constexpr size_t d = SpaceType::d;
  static const constexpr size_t r = SpaceType::r;
  static const constexpr size_t rC = SpaceType::rC;
  using R = typename SpaceType::R;

  DirichletConstraints(const BoundaryInfoType& bnd_info, const SpaceType& space, const bool set = true)
    : Propagator(this)
    , boundary_info_(bnd_info)
    , space_(space)
    , set_(set)
  {
  }

  void apply_local(const ElementType& element) override final
  {
    std::set<size_t> local_DoFs;
    const auto& fe = space_.finite_element(element.type());
    const auto& reference_element = ReferenceElements<double, d>::general(element.geometry().type());
    const auto local_key_indices = fe.coefficients().local_key_indices();
    const auto intersection_it_end = space_.grid_view().iend(element);
    for (auto intersection_it = space_.grid_view().ibegin(element); intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundary_info_.type(intersection) == XT::Grid::DirichletBoundary()
          || (!intersection.neighbor() && !intersection.boundary())) {
        const auto intersection_index = intersection.indexInInside();
        for (const auto& local_DoF : local_key_indices[1][intersection_index])
          local_DoFs.insert(local_DoF);
        for (int ii = 0; ii < reference_element.size(intersection_index, 1, d); ++ii) {
          const auto element_vertex_id = reference_element.subEntity(intersection_index, 1, ii, d);
          for (const auto& local_DoF : local_key_indices[d][element_vertex_id])
            local_DoFs.insert(local_DoF);
        }
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
        matrix.unit_col(DoF);
        vector[DoF] = 0.0;
      }
    } else {
      for (const auto& DoF : dirichlet_DoFs_) {
        matrix.clear_row(DoF);
        vector[DoF] = 0.0;
      }
    }
  } // ... apply(...)

  void finalize() override
  {
    this->finalize_imp();
  }

  BaseType* copy() override
  {
    return Propagator::copy_imp();
  }

  std::set<size_t> result() const
  {
    return dirichlet_DoFs_;
  }

  void set_result(std::set<size_t> res)
  {
    dirichlet_DoFs_ = res;
  }

private:
  const BoundaryInfoType& boundary_info_;
  const SpaceInterface<GridView, r, rC, R>& space_;
  const bool set_;
  std::set<size_t> dirichlet_DoFs_;
}; // class DirichletConstraints


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
