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

template <size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<HatFunctions<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = HatFunctions<double, 1, double, dimRange, 1, 1>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t quad_order = 15;
  static constexpr size_t num_quad_refinements = 0;

#if HAVE_CLP
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
struct RealizabilityLimiterChooser<PiecewiseMonomials<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = PiecewiseMonomials<double, 1, double, dimRange, 1, 1>;
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
struct RealizabilityLimiterChooser<RealSphericalHarmonics<double, double, order, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = RealSphericalHarmonics<double, double, order, 3>;
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
struct RealizabilityLimiterChooser<HatFunctions<double, 3, double, refinements, 1, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = HatFunctions<double, 3, double, refinements, 1, 3>;
  using LocalRealizabilityLimiterType =
      NonLimitingLocalRealizabilityLimiter<AnalyticalFluxType, DiscreteFunctionType, BasisfunctionType>;
  static constexpr size_t dimRange = BasisfunctionType::dimRange;
  static constexpr size_t quad_order = 7; // fekete rule number 7
  static constexpr size_t num_quad_refinements = 2;

#if HAVE_CLP
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
struct RealizabilityLimiterChooser<PiecewiseMonomials<double, 3, double, refinements, 1, 3, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using BasisfunctionType = PiecewiseMonomials<double, 3, double, refinements, 1, 3>;
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
struct SourceBeamPnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33066818456325014 : 0.33107004463413914;
  static constexpr double l2norm = reconstruct ? 0.4615751405564803 : 0.44609169128863851;
  static constexpr double linfnorm = reconstruct ? 1.1553979882432861 : 1.0882801946666156;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497237 : 0.33146794280839997;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337588411 : 0.33140398337567956;
  static constexpr double l2norm = reconstruct ? 0.4558335407458029 : 0.44484887611129575;
  static constexpr double linfnorm = reconstruct ? 0.99172184304968958 : 0.98930905293217597;
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
struct SourceBeamMnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140386483040757 : 0.33140386482516998;
  static constexpr double l2norm = reconstruct ? 0.45585375203639722 : 0.44485813651836886;
  static constexpr double linfnorm = reconstruct ? 0.99172157113121273 : 0.98930892899939982;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337608101 : 0.33140398337582355;
  static constexpr double l2norm = reconstruct ? 0.45585374774065163 : 0.44485696909271483;
  static constexpr double linfnorm = reconstruct ? 0.99172209692400415 : 0.98930944853186242;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337588411 : 0.33140398337567956;
  static constexpr double l2norm = reconstruct ? 0.4558335407458029 : 0.44484887611129575;
  static constexpr double linfnorm = reconstruct ? 0.99172184304968958 : 0.98930905293217597;
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
struct PlaneSourcePnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.9627559791618099 : 2.7793543802214402;
  static constexpr double linfnorm = reconstruct ? 7.5368337466833273 : 5.9468208917837284;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000557;
  static constexpr double l2norm = reconstruct ? 2.892587690555561 : 2.7677861047579322;
  static constexpr double linfnorm = reconstruct ? 6.9955083584307651 : 5.8898335510903852;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.000000024000026 : 2.0000000240000273;
  static constexpr double l2norm = reconstruct ? 2.881005248537496 : 2.7713504721240083;
  static constexpr double linfnorm = reconstruct ? 6.9331778582604997 : 6.0086435546642116;
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
struct PlaneSourceMnExpectedResults<LegendrePolynomials<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.7921993086492169 : 2.7461013585034388;
  static constexpr double linfnorm = reconstruct ? 4.8849177621513 : 5.3276983579096191;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<HatFunctions<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000158;
  static constexpr double l2norm = reconstruct ? 2.7968403961890758 : 2.7457197084995624;
  static constexpr double linfnorm = reconstruct ? 5.247368376105662 : 4.9918923122990027;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<PiecewiseMonomials<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000239999913 : 2.0000000239999904;
  static constexpr double l2norm = reconstruct ? 2.8215879031830924 : 2.7633864171093845;
  static constexpr double linfnorm = reconstruct ? 6.0674052799293623 : 6.2607864745536039;
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
struct PointSourcePnExpectedResults<RealSphericalHarmonics<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0007954640626406 : 1.0007954640534624;
  static constexpr double l2norm = reconstruct ? 2.6992308546885653 : 2.6970513662067956;
  static constexpr double linfnorm = reconstruct ? 10.379541276684977 : 10.414316177722712;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctions<double, 3, double, 0, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008094159849688 : 1.0008094159743741;
  static constexpr double l2norm = reconstruct ? 2.7092776186023921 : 2.7069983342698274;
  static constexpr double linfnorm = reconstruct ? 10.423991903881772 : 10.456911277964574;
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
struct PointSourcePnExpectedResults<HatFunctions<double, 3, double, 1, 1, 3>, reconstruct>
{
  // If Fekete is not available, we use a different quadrature, which gives slightly different results
  static_assert(!reconstruct, "Results with reconstruction not available yet!");
#if HAVE_FEKETE
  static constexpr double l1norm = 1.0007953665771843;
  static constexpr double l2norm = 2.7065005653281369;
  static constexpr double linfnorm = 10.456533271787738;
#else
  static constexpr double l1norm = 1.0008039111672102;
  static constexpr double l2norm = 2.7065211670792411;
  static constexpr double linfnorm = 10.456611570230423;
#endif
  // see above
  static constexpr double tol = reconstruct ? 1e-5 : 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PiecewiseMonomials<double, 3, double, 0, 1, 3>, reconstruct>
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
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PiecewiseMonomials<double, 3, double, 1, 1, 3>, reconstruct>
{
  static_assert(!reconstruct, "Results with reconstruction not available yet!");
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = 1.0007953665769933;
  static constexpr double l2norm = 2.7065732611564592;
  static constexpr double linfnorm = 10.457047161924061;
#else
  static constexpr double l1norm = 1.0008292531061092;
  static constexpr double l2norm = 2.7066524774407608;
  static constexpr double linfnorm = 10.457348661644719;
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

// PointSourceMn
template <class BasisfunctionImp, bool reconstruct>
struct PointSourceMnExpectedResults;

template <bool reconstruct>
struct PointSourceMnExpectedResults<RealSphericalHarmonics<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0007954640632573 : 1.0007954640632366;
  static constexpr double l2norm = reconstruct ? 2.6875983831354029 : 2.6817153143915298;
  static constexpr double linfnorm = reconstruct ? 10.360218210413363 : 10.365996649935104;
  static constexpr double tol = 1e-14;
};

template <bool reconstruct>
struct PointSourceMnExpectedResults<HatFunctions<double, 3, double, 0, 1, 3>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0007954632958449 : 1.0007954632958254;
  static constexpr double l2norm = reconstruct ? 2.6947689708516487 : 2.6892993875002693;
  static constexpr double linfnorm = reconstruct ? 10.379119469345591 : 10.395364963148149;
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
struct PointSourceMnExpectedResults<PiecewiseMonomials<double, 3, double, 0, 1, 3, 1>, reconstruct>
{
// If Fekete is not available, we use a different quadrature, which gives slightly different results
#if HAVE_FEKETE
  static constexpr double l1norm = reconstruct ? 1.0008094159885297 : 1.0008094159873213;
  static constexpr double l2norm = reconstruct ? 2.6984032921653527 : 2.6882435933317925;
  static constexpr double linfnorm = reconstruct ? 10.391298207973659 : 10.394254425939714;
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


} // namespace KineticTransport
} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH