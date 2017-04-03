// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_BASE_HH
#define DUNE_GDT_TEST_PROLONGATIONS_BASE_HH

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/la/container.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/expression.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/l2.hh>

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {

template <class GridType>
struct ProlongationOperatorsBaseGridHolder
{
  typedef Dune::XT::Grid::GridProvider<GridType> GridProviderType;

  ProlongationOperatorsBaseGridHolder()
    : grid_provider_(XT::Grid::make_cube_grid<GridType>(0.0, 1.0, 2u))
  {
    grid_provider_.global_refine(1);
  }

  GridProviderType grid_provider_;
}; // struct ProlongationOperatorsBaseGridHolder


template <class CoarseSpaceType, class FineSpaceType>
struct ProlongationOperatorsBase
    : public ::testing::Test,
      public ProlongationOperatorsBaseGridHolder<typename FineSpaceType::GridLayerType::Grid>
{
  typedef ProlongationOperatorsBaseGridHolder<typename FineSpaceType::GridLayerType::Grid> BaseType;
  typedef typename FineSpaceType::GridLayerType GridLayerType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename FineSpaceType::RangeFieldType RangeFieldType;
  static const size_t dimRange = FineSpaceType::dimRange;
  typedef Dune::XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef typename Dune::XT::LA::Container<RangeFieldType>::VectorType VectorType;
  typedef Dune::GDT::DiscreteFunction<CoarseSpaceType, VectorType> CoarseDiscreteFunctionType;
  typedef Dune::GDT::DiscreteFunction<FineSpaceType, VectorType> FineDiscreteFunctionType;

  static constexpr double default_tolerance = 1e-15;
  static constexpr double alugrid_tolerance = 3.8e-11;

  ProlongationOperatorsBase()
    : function_("x", "x[0]", 1, "function")
    , coarse_space_(grid_provider_.template level<CoarseSpaceType::part_view_type>(0))
    , fine_space_(grid_provider_.template level<FineSpaceType::part_view_type>(grid_provider_.grid().maxLevel()))
    , coarse_discrete_function_(coarse_space_, "coarse discrete function")
    , fine_discrete_function_(fine_space_, "fine discrete function")
    , prepared_(false)
  {
  }

  void prepare(const double tolerance)
  {
    if (prepared_)
      return;
    // first, project an anlytical function onto the coarse grid
    project(function_, coarse_discrete_function_);
    // since the projection operator is tested elsewhere we are confident this worked, but we check anyway
    const auto coarse_l2_error =
        make_l2_operator(coarse_space_.grid_layer(), 2)->induced_norm(function_ - coarse_discrete_function_);
    if (coarse_l2_error > tolerance)
      DUNE_THROW(Dune::XT::Common::Exceptions::internal_error,
                 "This should not happen, the L2 product and projection operators are tested elsewhere!\n"
                     << coarse_l2_error
                     << " vs. "
                     << tolerance);
    prepared_ = true;
  } // ... prepare(...)

  void produces_correct_results(const RangeFieldType& tolerance = default_tolerance)
  {
    prepare(tolerance);
    ProlongationOperatorType(fine_space_.grid_layer()).apply(coarse_discrete_function_, fine_discrete_function_);
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_layer(), 2)->induced_norm(function_ - fine_discrete_function_);
    EXPECT_LE(fine_l2_error, tolerance);
  } // ... produces_correct_results(...)

  using BaseType::grid_provider_;

  const FunctionType function_;
  const CoarseSpaceType coarse_space_;
  const FineSpaceType fine_space_;
  CoarseDiscreteFunctionType coarse_discrete_function_;
  FineDiscreteFunctionType fine_discrete_function_;
  bool prepared_;
}; // ProlongationOperatorsBase

template <class T, class U>
constexpr double ProlongationOperatorsBase<T, U>::default_tolerance;
template <class T, class U>
constexpr double ProlongationOperatorsBase<T, U>::alugrid_tolerance;

} // namespace internal


/**
 * \note Assumes that project and Products::L2 does the right thing.
 */
template <class CoarseSpaceType, class FineSpaceType, template <class, class, class> class ProlongationOperatorImp>
struct LocalizableProlongationOperatorBase : public internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType>
{
  typedef internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType> BaseType;
  using typename BaseType::GridLayerType;
  using typename BaseType::CoarseDiscreteFunctionType;
  using typename BaseType::FineDiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  typedef ProlongationOperatorImp<GridLayerType, CoarseDiscreteFunctionType, FineDiscreteFunctionType>
      ProlongationOperatorType;

  using BaseType::prepare;

  void produces_correct_results(const RangeFieldType& tolerance = double(BaseType::default_tolerance))
  {
    prepare(tolerance);
    ProlongationOperatorType(fine_space_.grid_layer(), coarse_discrete_function_, fine_discrete_function_).apply();
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_layer(), 2)->induced_norm(function_ - fine_discrete_function_);
    EXPECT_LE(fine_l2_error, tolerance);
  } // ... produces_correct_results(...)

  using BaseType::function_;
  using BaseType::fine_space_;
  using BaseType::coarse_discrete_function_;
  using BaseType::fine_discrete_function_;
}; // LocalizableProlongationOperatorBase


/**
 * \note Assumes that project and Products::L2 does the right thing.
 */
template <class CoarseSpaceType, class FineSpaceType, template <class, class> class ProlongationOperatorImp>
struct ProlongationOperatorBase : public internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType>
{
  typedef internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType> BaseType;
  using typename BaseType::GridLayerType;
  using typename BaseType::CoarseDiscreteFunctionType;
  using typename BaseType::FineDiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  typedef ProlongationOperatorImp<GridLayerType, RangeFieldType> ProlongationOperatorType;

  using BaseType::prepare;

  void produces_correct_results(const RangeFieldType& tolerance = double(BaseType::default_tolerance))
  {
    prepare(tolerance);
    ProlongationOperatorType op(fine_space_.grid_layer());
    op.apply(coarse_discrete_function_, fine_discrete_function_);
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_layer(), 2)->induced_norm(function_ - fine_discrete_function_);
    EXPECT_LE(fine_l2_error, tolerance);
  } // ... produces_correct_results(...)

  using BaseType::function_;
  using BaseType::fine_space_;
  using BaseType::coarse_discrete_function_;
  using BaseType::fine_discrete_function_;
}; // ProlongationOperatorBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_BASE_HH
