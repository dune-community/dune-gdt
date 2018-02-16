// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Tim Keil        (2017)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_ASSEMBLER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WRAPPER_HH

#include <type_traits>

#include <dune/xt/la/container/interfaces.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {
namespace internal {

// //////////////////////
// // wrap constraints //
// //////////////////////

template <class TestSpaceType, class AnsatzSpaceType, class GridLayerType, class ConstraintsType>
class ConstraintsWrapper : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(AlwaysFalse<ConstraintsType>::value, "Please add a specialization for these Constraints!");
};


// given DirichletConstraints

template <class TestSpaceType, class AnsatzSpaceType, class GridLayerType>
class ConstraintsWrapper<TestSpaceType,
                         AnsatzSpaceType,
                         GridLayerType,
                         DirichletConstraints<XT::Grid::extract_intersection_t<GridLayerType>>>
    : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "TestSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<AnsatzSpaceType>::value, "AnsatzSpaceType has to be derived from SpaceInterface!");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;
  typedef DirichletConstraints<XT::Grid::extract_intersection_t<GridLayerType>> ConstraintsType;

public:
  using typename BaseType::EntityType;

  ConstraintsWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                     const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                     const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                     ConstraintsType& constraints)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , constraints_(constraints)
    , thread_local_constraints_(constraints_.boundary_info(), constraints_.size())
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    test_space_->local_constraints(*ansatz_space_, entity, *thread_local_constraints_);
  }

  void finalize() override final
  {
    DUNE_UNUSED std::lock_guard<std::mutex> mutex_guard(constraints_.mutex_);
    constraints_.dirichlet_DoFs_.insert(thread_local_constraints_->dirichlet_DoFs_.begin(),
                                        thread_local_constraints_->dirichlet_DoFs_.end());
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  ConstraintsType& constraints_;
  XT::Common::PerThreadValue<ConstraintsType> thread_local_constraints_;
}; // class ConstraintsWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
