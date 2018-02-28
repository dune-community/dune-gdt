// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2013 - 2014, 2016 - 2018)
//   Tim Keil        (2017)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/assembler/functional-assemblers.hh>
#include <dune/gdt/local/assembler/two-form-assemblers.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class GridView,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class RangeField = double,
          class TestGridView = GridView,
          class AnsatzGridView = GridView,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols>
class GlobalAssembler : public XT::Grid::Walker<GridView>
{
  // no need to check the rest (is done in SpaceInterface)
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = GlobalAssembler<GridView,
                                   test_range_dim,
                                   test_range_dim_cols,
                                   RangeField,
                                   TestGridView,
                                   AnsatzGridView,
                                   ansatz_range_dim,
                                   ansatz_range_dim_cols>;
  using BaseType = XT::Grid::Walker<GridView>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;
  using TestSpaceType = SpaceInterface<TestGridView, test_range_dim, test_range_dim_cols, RangeField>;
  using AnsatzSpaceType = SpaceInterface<AnsatzGridView, ansatz_range_dim, ansatz_range_dim_cols, RangeField>;

  using ElementFilterType = XT::Grid::ElementFilter<GridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<GridViewType>;

  using LocalElementFunctionalType =
      LocalElementFunctionalInterface<ElementType, test_range_dim, test_range_dim_cols, RangeField>;

  using LocalElementTwoFormType = LocalElementTwoFormInterface<ElementType,
                                                               test_range_dim,
                                                               test_range_dim_cols,
                                                               RangeField,
                                                               RangeField,
                                                               ansatz_range_dim,
                                                               ansatz_range_dim_cols,
                                                               RangeField>;

  GlobalAssembler(GridViewType grd_vw, const TestSpaceType& test_sp, const AnsatzSpaceType& ansatz_sp)
    : BaseType(std::move(grd_vw))
    , test_space_(test_sp)
    , ansatz_space_(ansatz_sp)
  {
  }

  template <class GV,
            size_t r,
            size_t rC,
            class R, /* Only enable this ctor, iff */
            typename = typename std::enable_if</* the type of space is TestSpaceType and */ (
                                                   std::is_same<GV, TestGridView>::value && (r == test_range_dim)
                                                   && (rC == test_range_dim_cols)
                                                   && std::is_same<R, RangeField>::value)
                                               && /* if ansatz and test space type coincide. */ std::
                                                      is_same<AnsatzSpaceType, TestSpaceType>::value>::type>
  GlobalAssembler(GridViewType grd_vw, const SpaceInterface<GV, r, rC, R>& space)
    : GlobalAssembler(std::move(grd_vw), space, space)
  {
  }

  template <class GV,
            size_t r,
            size_t rC,
            class R, /* Only enable this ctor, iff */
            typename = typename std::enable_if</* the type of space is TestSpaceType, */ (
                                                   std::is_same<GV, TestGridView>::value && (r == test_range_dim)
                                                   && (rC == test_range_dim_cols)
                                                   && std::is_same<R, RangeField>::value)
                                               && /* ansatz and test space type coincide and */ std::
                                                      is_same<AnsatzSpaceType, TestSpaceType>::value
                                               && /* the grid view type of space and GridViewType coincide */
                                               std::is_same<GV, GridViewType>::value>::type>
  GlobalAssembler(const SpaceInterface<GV, r, rC, R>& space)
    : GlobalAssembler(space.grid_view(), space, space)
  {
  }

  GlobalAssembler(const ThisType& other) = delete;
  GlobalAssembler(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  using BaseType::append;

  template <class V>
  ThisType& append(const LocalElementFunctionalType& local_functional,
                   XT::LA::VectorInterface<V>& global_vector,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    using LocalAssemblerType = LocalElementFunctionalAssembler<typename XT::LA::VectorInterface<V>::derived_type,
                                                               GridViewType,
                                                               test_range_dim,
                                                               test_range_dim_cols,
                                                               RangeField,
                                                               TestGridView>;
    this->append(new LocalAssemblerType(test_space_, local_functional, global_vector.as_imp(), param), filter);
    return *this;
  }

  template <class M>
  ThisType& append(const LocalElementTwoFormType& local_two_form,
                   XT::LA::MatrixInterface<M>& global_matrix,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    using LocalAssemblerType = LocalElementTwoFormAssembler<typename XT::LA::MatrixInterface<M>::derived_type,
                                                            GridViewType,
                                                            test_range_dim,
                                                            test_range_dim_cols,
                                                            RangeField,
                                                            TestGridView,
                                                            AnsatzGridView,
                                                            ansatz_range_dim,
                                                            ansatz_range_dim_cols>;
    this->append(new LocalAssemblerType(test_space_, ansatz_space_, local_two_form, global_matrix.as_imp(), param),
                 filter);
    return *this;
  }

  void assemble(const bool use_tbb = false, const bool clear_functors = true)
  {
    this->walk(use_tbb, clear_functors);
  }

  template <class Partitioning>
  void assemble(Partitioning& partitioning, const bool clear_functors = true)
  {
    this->walk(partitioning, clear_functors);
  }

protected:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
}; // class GlobalAssembler


template <class GridView,
          class TestGridView,
          size_t t_r,
          size_t t_rC,
          class R,
          class AnsatzGridView,
          size_t a_r,
          size_t a_rC>
GlobalAssembler<GridView, t_r, t_rC, R, TestGridView, AnsatzGridView, a_r, a_rC>
make_global_assembler(GridView grid_view,
                      const SpaceInterface<TestGridView, t_r, t_rC, R>& test_space,
                      const SpaceInterface<AnsatzGridView, a_r, a_rC, R>& ansatz_space,
                      const SpaceInterface<TestGridView, t_r, t_rC, R>& outer_test_space,
                      const SpaceInterface<AnsatzGridView, a_r, a_rC, R>& outer_ansatz_space)
{
  return GlobalAssembler<GridView, t_r, t_rC, R, TestGridView, AnsatzGridView, a_r, a_rC>(
      grid_view, test_space, ansatz_space, outer_test_space, outer_ansatz_space);
}


template <class GridView,
          class TestGridView,
          size_t t_r,
          size_t t_rC,
          class R,
          class AnsatzGridView,
          size_t a_r,
          size_t a_rC>
GlobalAssembler<GridView, t_r, t_rC, R, TestGridView, AnsatzGridView, a_r, a_rC>
make_global_assembler(GridView grid_view,
                      const SpaceInterface<TestGridView, t_r, t_rC, R>& test_space,
                      const SpaceInterface<AnsatzGridView, a_r, a_rC, R>& ansatz_space)
{
  return GlobalAssembler<GridView, t_r, t_rC, R, TestGridView, AnsatzGridView, a_r, a_rC>(
      grid_view, test_space, ansatz_space);
}


template <class GridView, class SpaceGridView, size_t r, size_t rC, class R>
GlobalAssembler<GridView, r, rC, R, SpaceGridView>
make_global_assembler(GridView grid_view, const SpaceInterface<SpaceGridView, r, rC, R>& space)
{
  return GlobalAssembler<GridView, r, rC, R, SpaceGridView>(grid_view, space);
}


template <class GV, size_t r, size_t rC, class R>
GlobalAssembler<GV, r, rC, R> make_global_assembler(const SpaceInterface<GV, r, rC, R>& space)
{
  return GlobalAssembler<GV, r, rC, R>(space);
}


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
