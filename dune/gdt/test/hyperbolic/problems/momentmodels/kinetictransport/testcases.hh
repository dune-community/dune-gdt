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
#include <dune/gdt/operators/fv/entropybased/realizability.hh>
#include <dune/gdt/operators/fv/reconstruction/slopes.hh>

#include "planesource.hh"
#include "pointsource.hh"
#include "sourcebeam.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {
namespace KineticTransport {


// choose RealizabilityLimiter suitable for BasisfunctionImp
template <class BasisfunctionImp, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser;

#if HAVE_CLP
template <size_t order, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<LegendrePolynomials<double, double, order>, AnalyticalFluxType, DiscreteFunctionType>
{
  using BasisfunctionType = LegendrePolynomials<double, double, order>;
  //  using LocalRealizabilityLimiterType =
  //      ClpLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;

  template <class MatrixType, class QuadratureType>
  static std::unique_ptr<LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>>
  make_slope(const BasisfunctionType& basis_functions, const QuadratureType& quadrature, const double epsilon)
  {
    using SlopeType = LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>;
    return std::make_unique<SlopeType>(basis_functions, quadrature, epsilon);
  }
};
#endif

template <size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<HatFunctions<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = HatFunctions<double, 1, double, dimRange, 1, 1>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;

#if HAVE_CLP
  template <class MatrixType, class QuadratureType>
  static std::unique_ptr<LpPositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const QuadratureType& /*quadrature*/, const double epsilon)
  {
    using SlopeType = LpPositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#else // HAVE_CLP
  template <class MatrixType, class QuadratureType>
  static std::unique_ptr<PositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const QuadratureType& /*quadrature*/, const double epsilon)
  {
    using SlopeType = PositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#endif // HAVE_CLP
};

template <size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<PiecewiseMonomials<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = PiecewiseMonomials<double, 1, double, dimRange, 1, 1>;
  using LocalRealizabilityLimiterType =
      DgLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;

  template <class MatrixType, class QuadratureType>
  static std::unique_ptr<Dg1dRealizabilityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const QuadratureType& /*quadrature*/, const double epsilon)
  {
    using SlopeType = Dg1dRealizabilityLimitedSlope<double, dimRange, MatrixType>;
    //   using SlopeType = MinmodSlope<RangeType, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
};


// SourceBeam Pn
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


// SourceBeam Mn
template <class BasisfunctionImp, bool reconstruct>
struct SourceBeamMnExpectedResults;

template <bool reconstruct>
struct SourceBeamMnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33066818456325014 : 0.33107004463413914;
  static constexpr double l2norm = reconstruct ? 0.4615751405564803 : 0.44609169128863851;
  static constexpr double linfnorm = reconstruct ? 1.1553979882432861 : 1.0882801946666156;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497237 : 0.33146794280839997;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337610927 : 0.33140398337603194;
  static constexpr double l2norm = reconstruct ? 0.47294828933204158 : 0.45667075585121392;
  static constexpr double linfnorm = reconstruct ? 1.0490804598503622 : 0.99004736850989217;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct SourceBeamMnTestCase
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
  using ProblemType = SourceBeamMn<BasisfunctionType, GridLayerType, DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = SourceBeamMnExpectedResults<BasisfunctionImp, reconstruction>;
  using RealizabilityLimiterChooserType =
      RealizabilityLimiterChooser<BasisfunctionType, typename ProblemType::FluxType, DiscreteFunctionType>;
};


// PlaneSource Pn
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


// PlaneSource Mn
template <class BasisfunctionImp, bool reconstruct>
struct PlaneSourceMnExpectedResults;

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.9627559791618099 : 2.7793543802214402;
  static constexpr double linfnorm = reconstruct ? 7.5368337466833273 : 5.9468208917837284;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000557;
  static constexpr double l2norm = reconstruct ? 2.892587690555561 : 2.7677861047579322;
  static constexpr double linfnorm = reconstruct ? 6.9955083584307651 : 5.8898335510903852;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.000000024000026 : 2.0000000240000273;
  static constexpr double l2norm = reconstruct ? 2.881005248537496 : 2.7713504721240083;
  static constexpr double linfnorm = reconstruct ? 6.9331778582604997 : 6.0086435546642116;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct PlaneSourceMnTestCase : SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using typename BaseType::DiscreteFunctionType;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using ProblemType = PlaneSourceMn<BasisfunctionImp, typename BaseType::GridLayerType, DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PlaneSourceMnExpectedResults<BasisfunctionImp, reconstruction>;
  using RealizabilityLimiterChooserType =
      RealizabilityLimiterChooser<BasisfunctionImp, typename ProblemType::FluxType, DiscreteFunctionType>;
};


// PointSourcePn
template <class BasisfunctionImp, bool reconstruct>
struct PointSourcePnExpectedResults;

template <bool reconstruct>
struct PointSourcePnExpectedResults<RealSphericalHarmonics<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0007954640626406 : 1.0007954640534238;
  static constexpr double l2norm = reconstruct ? 2.7177565161122006 : 2.7163083579825025;
  static constexpr double linfnorm = reconstruct ? 10.461558474249745 : 10.498572083981468;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctions<double, 3, double, 6, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008094159849688 : 1.0008094159743741;
  static constexpr double l2norm = reconstruct ? 2.7092776186023921 : 2.7069983342698274;
  static constexpr double linfnorm = reconstruct ? 10.423991903881772 : 10.456911277964574;
#else
  static constexpr double l1norm = reconstruct ? 1.0008292531174403 : 1.0008292531057066;
  static constexpr double l2norm = reconstruct ? 2.7095626432312425 : 2.7070581236565103;
  static constexpr double linfnorm = reconstruct ? 10.424226802303412 : 10.457145890791487;
#endif
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PiecewiseMonomials<double, 3, double, 32, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008094159850585 : 1.000809415974838;
  static constexpr double l2norm = reconstruct ? 2.7098602740535496 : 2.7065939033692201;
  static constexpr double linfnorm = reconstruct ? 10.427604575554344 : 10.457121881221033;
#else
  static constexpr double l1norm = reconstruct ? 1.0008292531175822 : 1.0008292531061092;
  static constexpr double l2norm = reconstruct ? 2.7099187578817849 : 2.7066524774407608;
  static constexpr double linfnorm = reconstruct ? 10.427830136315574 : 10.457348661644719;
#endif
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
