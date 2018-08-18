// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/gdt/timestepper/interface.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

#include "sourcebeam.hh"
#include "planesource.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


template <class BasisfunctionImp, bool reconstruct>
struct SourceBeamPnExpectedResults;

template <bool reconstruct>
struct SourceBeamPnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497237 : 0.33146794280839997;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct SourceBeamPnTestCase
{
  using BasisfunctionType = BasisfunctionImp;
  static constexpr size_t dimDomain = BasisfunctionType::dimDomain;
  static constexpr size_t dimRange = BasisfunctionType::dimRange;
  static constexpr auto time_stepper_method = TimeStepperMethods::explicit_rungekutta_second_order_ssp;
  static constexpr auto rhs_time_stepper_method = TimeStepperMethods::matrix_exponential;
  static constexpr auto time_stepper_splitting_method = TimeStepperSplittingMethods::strang;
  using DomainFieldType = typename BasisfunctionType::DomainFieldType;
  using RangeFieldType = typename BasisfunctionType::RangeFieldType;
  using GridType = GridImp;
  using GridLayerType = typename GridType::LeafGridView;
  using SpaceType = FvProductSpace<GridLayerType, RangeFieldType, dimRange, 1>;
  using VectorType = typename Dune::XT::LA::Container<RangeFieldType, Dune::XT::LA::default_backend>::VectorType;
  using DiscreteFunctionType = DiscreteFunction<SpaceType, VectorType>;
  using ProblemType = SourceBeamPn<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = SourceBeamPnExpectedResults<BasisfunctionImp, reconstruction>;
};

template <class BasisfunctionImp, bool reconstruct>
struct PlaneSourcePnExpectedResults;

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000557;
  static constexpr double l2norm = reconstruct ? 2.892587690555561 : 2.7677861047579322;
  static constexpr double linfnorm = reconstruct ? 6.9955083584307651 : 5.8898335510903852;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct PlaneSourcePnTestCase : SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using ProblemType =
      PlaneSourcePn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PlaneSourcePnExpectedResults<BasisfunctionImp, reconstruction>;
};


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH
