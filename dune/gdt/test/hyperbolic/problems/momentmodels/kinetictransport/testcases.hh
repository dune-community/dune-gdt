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

#include "planesource.hh"
#include "pointsource.hh"
#include "sourcebeam.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


// SourceBeam
template <class BasisfunctionImp, bool reconstruct>
struct SourceBeamPnExpectedResults;

template <bool reconstruct>
struct SourceBeamPnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33066818456325014 : 0.33107004463413914;
  static constexpr double l2norm = reconstruct ? 0.4615751405564803 : 0.44609169128863851;
  static constexpr double linfnorm = reconstruct ? 1.1553979882432861 : 1.0882801946666156;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497237 : 0.33146794280839997;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337610927 : 0.33140398337603194;
  static constexpr double l2norm = reconstruct ? 0.47294828933204158 : 0.45667075585121392;
  static constexpr double linfnorm = reconstruct ? 1.0490804598503622 : 0.99004736850989217;
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

// PlaneSource
template <class BasisfunctionImp, bool reconstruct>
struct PlaneSourcePnExpectedResults;

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.9627559791618099 : 2.7793543802214402;
  static constexpr double linfnorm = reconstruct ? 7.5368337466833273 : 5.9468208917837284;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000557;
  static constexpr double l2norm = reconstruct ? 2.892587690555561 : 2.7677861047579322;
  static constexpr double linfnorm = reconstruct ? 6.9955083584307651 : 5.8898335510903852;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.000000024000026 : 2.0000000240000273;
  static constexpr double l2norm = reconstruct ? 2.881005248537496 : 2.7713504721240083;
  static constexpr double linfnorm = reconstruct ? 6.9331778582604997 : 6.0086435546642116;
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


// PointSource
template <class BasisfunctionImp, bool reconstruct>
struct PointSourcePnExpectedResults;

// template <bool reconstruct>
// struct PointSourcePnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
//{
//  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
//  static constexpr double l2norm = reconstruct ? 2.9627559791618099 : 2.7793543802214402;
//  static constexpr double linfnorm = reconstruct ? 7.5368337466833273 : 5.9468208917837284;
//};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctions<double, 3, double, 6, 1, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0029611747117857 : 1.0029611746535056;
  static constexpr double l2norm = reconstruct ? 2.714746717141042 : 2.7122234140989883;
  static constexpr double linfnorm = reconstruct ? 10.443809507911537 : 10.476793054121075;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PiecewiseMonomials<double, 3, double, 32, 1, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0088250566247499 : 1.0029611746514546;
  static constexpr double l2norm = reconstruct ? 2.7626619532892951 : 2.7118027445930162;
  static constexpr double linfnorm = reconstruct ? 10.787213242460366 : 10.476921363773437;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct PointSourcePnTestCase : SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using ProblemType =
      PointSourcePn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourcePnExpectedResults<BasisfunctionImp, reconstruction>;
};


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH
