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

#include "checkerboard.hh"
#include "planesource.hh"
#include "pointsource.hh"
#include "shadow.hh"
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
struct RealizabilityLimiterChooser<LegendreMomentBasis<double, double, order>, AnalyticalFluxType, DiscreteFunctionType>
{
  using BasisfunctionType = LegendreMomentBasis<double, double, order>;
  //  using LocalRealizabilityLimiterType =
  //      ClpLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 31;
  static constexpr size_t num_quad_refinements = 6;

  template <class MatrixType>
  static std::unique_ptr<LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>>
  make_slope(const BasisfunctionType& basis_functions, const double epsilon)
  {
    using SlopeType = LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>;
    return std::make_unique<SlopeType>(basis_functions, epsilon);
  }

  static std::unique_ptr<BasisfunctionType> make_basis_functions()
  {
    return std::make_unique<BasisfunctionType>(31, 6);
  }
};
#endif

#ifndef USE_LP_POSITIVITY_LIMITER
#define USE_LP_POSITIVITY_LIMITER 0
#endif // USE_LP_POSITIVITY_LIMITER
template <size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<HatFunctionMomentBasis<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = HatFunctionMomentBasis<double, 1, double, dimRange, 1, 1>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 15;
  static constexpr size_t num_quad_refinements = 0;

#if HAVE_CLP && USE_LP_POSITIVITY_LIMITER
  template <class MatrixType>
  static std::unique_ptr<LpPositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = LpPositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#else // HAVE_CLP
  template <class MatrixType>
  static std::unique_ptr<PositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = PositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#endif // HAVE_CLP
};

template <size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<PartialMomentBasis<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = PartialMomentBasis<double, 1, double, dimRange, 1, 1>;
  //  using LocalRealizabilityLimiterType =
  //      DgLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 15;
  static constexpr size_t num_quad_refinements = 0;

  template <class MatrixType>
  static std::unique_ptr<Dg1dRealizabilityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& basis_functions, const double epsilon)
  {
    using SlopeType = Dg1dRealizabilityLimitedSlope<double, dimRange, MatrixType>;
    //   using SlopeType = MinmodSlope<RangeType, MatrixType>;
    return std::make_unique<SlopeType>(basis_functions, epsilon);
  }

  static std::unique_ptr<BasisfunctionType> make_basis_functions()
  {
    return std::make_unique<BasisfunctionType>();
  }
};

#if HAVE_CLP
template <size_t order, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<RealSphericalHarmonicsMomentBasis<double, double, order, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = RealSphericalHarmonicsMomentBasis<double, double, order, 3>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 2 * order + 2; // fekete rule number 7
  static constexpr size_t num_quad_refinements = 0;

  template <class MatrixType>
  static std::unique_ptr<LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>>
  make_slope(const BasisfunctionType& basis_functions, const double epsilon)
  {
    using SlopeType = LpConvexhullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>;
    return std::make_unique<SlopeType>(basis_functions, epsilon);
  }

  static std::unique_ptr<BasisfunctionType> make_basis_functions()
  {
    return std::make_unique<BasisfunctionType>();
  }
};
#endif

template <size_t refinements, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<HatFunctionMomentBasis<double, 3, double, refinements, 1, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = HatFunctionMomentBasis<double, 3, double, refinements, 1, 3>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t dimRange = BasisfunctionType::dimRange;
  static constexpr size_t quad_order = 7; // fekete rule number 7
  static constexpr size_t num_quad_refinements = 2;

#if HAVE_CLP && USE_LP_POSITIVITY_LIMITER
  template <class MatrixType>
  static std::unique_ptr<LpPositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = LpPositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#else // HAVE_CLP
  template <class MatrixType>
  static std::unique_ptr<PositivityLimitedSlope<double, dimRange, MatrixType>>
  make_slope(const BasisfunctionType& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = PositivityLimitedSlope<double, dimRange, MatrixType>;
    return std::make_unique<SlopeType>(epsilon);
  }
#endif // HAVE_CLP
};

#if HAVE_QHULL
template <size_t refinements, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<PartialMomentBasis<double, 3, double, refinements, 1, 3, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = PartialMomentBasis<double, 3, double, refinements, 1, 3>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 3; // fekete rule number 3
  static constexpr size_t num_quad_refinements = 0;

  template <class MatrixType>
  static std::unique_ptr<DgConvexHullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>>
  make_slope(const BasisfunctionType& basis_functions, const double epsilon)
  {
    using SlopeType = DgConvexHullRealizabilityLimitedSlope<BasisfunctionType, MatrixType>;
    //   using SlopeType = MinmodSlope<RangeType, MatrixType>;
    return std::make_unique<SlopeType>(basis_functions, epsilon);
  }
};
#endif

// SourceBeam Pn
template <class BasisfunctionImp, bool reconstruct>
struct SourceBeamPnExpectedResults;

template <bool reconstruct>
struct SourceBeamPnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33066818456325631 : 0.33107004463414536;
  static constexpr double l2norm = reconstruct ? 0.46157514055648519 : 0.44609169128864312;
  static constexpr double linfnorm = reconstruct ? 1.1553979882432861 : 1.0882801946666156;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497681 : 0.33146794280840425;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337610927 : 0.33140398337603194;
  static constexpr double l2norm = reconstruct ? 0.47294828933204164 : 0.45667075585121392;
  static constexpr double linfnorm = reconstruct ? 1.0490804598503625 : 0.99004736850989217;
  static constexpr double tol = 1e-14;
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
struct SourceBeamMnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140386482977746 : 0.33140386482489148;
  static constexpr double l2norm = reconstruct ? 0.45585374382359678 : 0.44485813652094131;
  static constexpr double linfnorm = reconstruct ? 0.99172157116705117 : 0.98930892901211287;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337583815 : 0.33140398337577592;
  static constexpr double l2norm = reconstruct ? 0.45585374774346554 : 0.44485696909274419;
  static constexpr double linfnorm = reconstruct ? 0.99172209693226199 : 0.98930944853353908;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.3314039833756896 : 0.33140398337567956;
  static constexpr double l2norm = reconstruct ? 0.45583354074698401 : 0.44484887611137452;
  static constexpr double linfnorm = reconstruct ? 0.99172184304767108 : 0.98930905293217597;
  static constexpr double tol = 1e-14;
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
struct PlaneSourcePnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.9616518419466558 : 2.7792352623482848;
  static constexpr double linfnorm = reconstruct ? 7.5355813391308644 : 5.9472849007944166;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000149;
  static constexpr double l2norm = reconstruct ? 2.8915349919892397 : 2.7676677008555917;
  static constexpr double linfnorm = reconstruct ? 6.9950740716997668 : 5.8904604670932663;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000239999896 : 2.0000000239999918;
  static constexpr double l2norm = reconstruct ? 2.8799152602279068 : 2.771228836660768;
  static constexpr double linfnorm = reconstruct ? 6.9320887958307775 : 6.0090382693364512;
  static constexpr double tol = 1e-14;
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
struct PlaneSourceMnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.7919962607432942 : 2.7461013585128331;
  static constexpr double linfnorm = reconstruct ? 4.8942859838731865 : 5.3276983579131185;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000158;
  static constexpr double l2norm = reconstruct ? 2.7968403961895407 : 2.745719708499653;
  static constexpr double linfnorm = reconstruct ? 5.2473683761050234 : 4.9918923122986048;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000239999913 : 2.0000000239999904;
  static constexpr double l2norm = reconstruct ? 2.8215879031862658 : 2.7633864171098814;
  static constexpr double linfnorm = reconstruct ? 6.0674052799283675 : 6.2607864745505113;
  static constexpr double tol = 1e-14;
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
struct PointSourcePnExpectedResults
{
  static constexpr double l1norm = reconstruct ? 1.0007954640626406 : 1.0007954640534238;
  static constexpr double l2norm = reconstruct ? 2.7177565161122006 : 2.7163083579825025;
  static constexpr double linfnorm = reconstruct ? 10.461558474249745 : 10.498572083981468;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0007954640626138 : 1.0007954640529835;
  static constexpr double l2norm = reconstruct ? 2.7018734377354714 : 2.6999617861637795;
  static constexpr double linfnorm = reconstruct ? 10.391017138477658 : 10.426558946481034;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008081476461759 : 1.0008081476360002;
  static constexpr double l2norm = reconstruct ? 2.7094988565517113 : 2.7073931584332405;
  static constexpr double linfnorm = reconstruct ? 10.424578957824737 : 10.458114243519278;
#else
  static constexpr double l1norm = reconstruct ? 1.0008292531174403 : 1.0008261813654091;
  static constexpr double l2norm = reconstruct ? 2.7097574564184925 : 2.7074480471822433;
  static constexpr double linfnorm = reconstruct ? 10.424798568201023 : 10.458329938780505;
#endif
  // The matrices in this test case all have eigenvalues [+-0.808311035811965, 0, 0, 0, 0].
  // Thus, the eigenvectors are not unique, and the eigensolvers are extremely sensitive
  // to numerical errors. A difference of 1e-16 in the jacobians entries suffices to
  // result in completely different eigenvectors. In all cases, the eigenvectors are
  // valid eigenvectors to the correct eigenvalues. However, this difference in the
  // eigendecomposition leads to differences in the results with linear reconstruction
  // that are larger than would be expected by pure numerical errors.
  static constexpr double tol = reconstruct ? 1e-5 : 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctionMomentBasis<double, 3, double, 1, 1, 3>, reconstruct>
{
  // If Fekete is not available, we use a different quadrature, which gives slightly different results
  static_assert(!reconstruct, "Results with reconstruction not available yet!");
#if HAVE_FEKETE
  static constexpr double l1norm = 1.0007953754379604;
  static constexpr double l2norm = 2.7069180208261652;
  static constexpr double linfnorm = 10.4578353362239;
#else
  static constexpr double l1norm = 1.0008031432622411;
  static constexpr double l2norm = 2.7069365027580932;
  static constexpr double linfnorm = 10.457905431746271;
#endif
  // see above
  static constexpr double tol = reconstruct ? 1e-5 : 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PartialMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008081476462325 : 1.0008081476365056;
  static constexpr double l2norm = reconstruct ? 2.7100893563191693 : 2.7070071430786014;
  static constexpr double linfnorm = reconstruct ? 10.428361386931211 : 10.458405903819866;
#else
  static constexpr double l1norm = reconstruct ? 1.0008261813756092 : 1.0008261813658712;
  static constexpr double l2norm = reconstruct ? 2.7101429263071872 : 2.7070608309129254;
  static constexpr double linfnorm = reconstruct ? 10.428568071618781 : 10.45861384878939;
#endif
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PartialMomentBasis<double, 3, double, 1, 1, 3>, reconstruct>
{
  static_assert(!reconstruct, "Results with reconstruction not available yet!");
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = 1.0007953754377614;
  static constexpr double l2norm = 2.7069901509556513;
  static constexpr double linfnorm = 10.458345325407425;
#else
  static constexpr double l1norm = 1.0008031432620421;
  static constexpr double l2norm = 2.707008617216633;
  static constexpr double linfnorm = 10.458415332332239;
#endif
  static constexpr double tol = 1e-14;
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

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct CheckerboardPnTestCase : SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using ProblemType =
      CheckerboardPn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourcePnExpectedResults<BasisfunctionImp, reconstruction>;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct ShadowPnTestCase : SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using ProblemType =
      ShadowPn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourcePnExpectedResults<BasisfunctionImp, reconstruction>;
};


// PointSourceMn
template <class BasisfunctionImp, bool reconstruct>
struct PointSourceMnExpectedResults;

template <bool reconstruct>
struct PointSourceMnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.000795464063172 : 1.0007954640599357;
  static constexpr double l2norm = reconstruct ? 2.6987974255774256 : 2.7004427948819156;
  static constexpr double linfnorm = reconstruct ? 10.392787269614516 : 10.429333233979499;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourceMnExpectedResults<HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0007954632958449 : 1.0007954633622282;
  static constexpr double l2norm = reconstruct ? 2.7073123070405787 : 2.7080473048720122;
  static constexpr double linfnorm = reconstruct ? 10.420529174853563 : 10.459179148064445;
#else
  static constexpr double l1norm = reconstruct ? 1.0008292531174403 : 1.0008292531057066;
  static constexpr double l2norm = reconstruct ? 2.7095647696183893 : 2.7070581236565103;
  static constexpr double linfnorm = reconstruct ? 10.424228258642177 : 10.457145890791487;
#endif
  // The matrices in this test case all have eigenvalues [+-0.808311035811965, 0, 0, 0, 0].
  // Thus, the eigenvectors are not unique, and the eigensolvers are extremely sensitive
  // to numerical errors. A difference of 1e-16 in the jacobians entries suffices to
  // result in completely different eigenvectors. In all cases, the eigenvectors are
  // valid eigenvectors to the correct eigenvalues. However, this difference in the
  // eigendecomposition leads to differences in the results with linear reconstruction
  // that are larger than would be expected by pure numerical errors.
  static constexpr double tol = reconstruct ? 1e-5 : 1e-14;
};

template <bool reconstruct>
struct PointSourceMnExpectedResults<PartialMomentBasis<double, 3, double, 0, 1, 3, 1>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.00080814764654 : 1.0008081476364998;
  static constexpr double l2norm = reconstruct ? 2.7100989741978259 : 2.7070280881780477;
  static constexpr double linfnorm = reconstruct ? 10.428410325022771 : 10.458452422334867;
#else
  static constexpr double l1norm = reconstruct ? 1.0008292531175822 : 1.0008292531061092;
  static constexpr double l2norm = reconstruct ? 2.7099187578817849 : 2.7066524774407608;
  static constexpr double linfnorm = reconstruct ? 10.427830136315574 : 10.457348661644719;
#endif
  static constexpr double tol = 1e-14;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct PointSourceMnTestCase : SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using ProblemType =
      PointSourceMn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourceMnExpectedResults<BasisfunctionImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<BasisfunctionImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct CheckerboardMnTestCase : SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using ProblemType =
      CheckerboardMn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourceMnExpectedResults<BasisfunctionImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<BasisfunctionImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};

template <class GridImp, class BasisfunctionImp, bool reconstruct>
struct ShadowMnTestCase : SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, BasisfunctionImp, reconstruct>;
  using ProblemType =
      ShadowMn<BasisfunctionImp, typename BaseType::GridLayerType, typename BaseType::DiscreteFunctionType>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourceMnExpectedResults<BasisfunctionImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<BasisfunctionImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH
