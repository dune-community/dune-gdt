// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_BASE_HH
#define DUNE_GDT_TEST_PROLONGATIONS_BASE_HH

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/test/gtest/gtest.h>

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
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;

  ProlongationOperatorsBaseGridHolder()
    : grid_provider_(0.0, 1.0, 2u)
  {
    grid_provider_.global_refine(1);
  }

  GridProviderType grid_provider_;
}; // struct ProlongationOperatorsBaseGridHolder


template <class CoarseSpaceType, class FineSpaceType>
struct ProlongationOperatorsBase
    : public ::testing::Test,
      public ProlongationOperatorsBaseGridHolder<typename FineSpaceType::GridViewType::Grid>
{
  typedef ProlongationOperatorsBaseGridHolder<typename FineSpaceType::GridViewType::Grid> BaseType;
  typedef typename FineSpaceType::GridViewType GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename FineSpaceType::RangeFieldType RangeFieldType;
  static const size_t dimRange = FineSpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef typename Dune::Stuff::LA::Container<RangeFieldType>::VectorType VectorType;
  typedef Dune::GDT::DiscreteFunction<CoarseSpaceType, VectorType> CoarseDiscreteFunctionType;
  typedef Dune::GDT::DiscreteFunction<FineSpaceType, VectorType> FineDiscreteFunctionType;

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
        make_l2_operator(coarse_space_.grid_view(), 2)->induced_norm(function_ - coarse_discrete_function_);
    if (coarse_l2_error > tolerance)
      DUNE_THROW(Dune::Stuff::Exceptions::internal_error,
                 "This should not happen, the L2 product and projection operators are tested elsewhere!\n"
                     << coarse_l2_error
                     << " vs. "
                     << tolerance);
    prepared_ = true;
  } // ... prepare(...)

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    prepare(tolerance);
    ProlongationOperatorType(fine_space_.grid_view()).apply(coarse_discrete_function_, fine_discrete_function_);
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_view(), 2)->induced_norm(function_ - fine_discrete_function_);
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


} // namespace internal


/**
 * \note Assumes that project and Products::L2 does the right thing.
 */
template <class CoarseSpaceType, class FineSpaceType, template <class, class, class> class ProlongationOperatorImp>
struct LocalizableProlongationOperatorBase : public internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType>
{
  typedef internal::ProlongationOperatorsBase<CoarseSpaceType, FineSpaceType> BaseType;
  using typename BaseType::GridViewType;
  using typename BaseType::CoarseDiscreteFunctionType;
  using typename BaseType::FineDiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  typedef ProlongationOperatorImp<GridViewType, CoarseDiscreteFunctionType, FineDiscreteFunctionType>
      ProlongationOperatorType;

  using BaseType::prepare;

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    prepare(tolerance);
    ProlongationOperatorType(fine_space_.grid_view(), coarse_discrete_function_, fine_discrete_function_).apply();
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_view(), 2)->induced_norm(function_ - fine_discrete_function_);
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
  using typename BaseType::GridViewType;
  using typename BaseType::CoarseDiscreteFunctionType;
  using typename BaseType::FineDiscreteFunctionType;
  using typename BaseType::RangeFieldType;
  typedef ProlongationOperatorImp<GridViewType, RangeFieldType> ProlongationOperatorType;

  using BaseType::prepare;

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    prepare(tolerance);
    ProlongationOperatorType op(fine_space_.grid_view());
    op.apply(coarse_discrete_function_, fine_discrete_function_);
    const auto fine_l2_error =
        make_l2_operator(fine_space_.grid_view(), 2)->induced_norm(function_ - fine_discrete_function_);
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
