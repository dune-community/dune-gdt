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

#include <dune/gdt/momentmodels/basisfunctions.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/tools/timestepper/interface.hh>
#include <dune/gdt/operators/reconstruction/slopes.hh>

#include "checkerboard.hh"
#include "planesource.hh"
#include "pointsource.hh"
#include "shadow.hh"
#include "sourcebeam.hh"

namespace Dune {
namespace GDT {


// choose RealizabilityLimiter suitable for MomentBasisImp
template <class GV, class MomentBasisImp, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser;

#if HAVE_CLP
template <class GV, size_t order, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   LegendreMomentBasis<double, double, order>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = LegendreMomentBasis<double, double, order>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t quad_order = 31;
  static constexpr size_t num_quad_refinements = 6;

  template <class EigenVectorWrapperType>
  static std::unique_ptr<LpConvexhullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& basis_functions, const double epsilon)
  {
    using SlopeType = LpConvexhullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, basis_functions, epsilon);
  }
};
#endif

#ifndef USE_LP_POSITIVITY_LIMITER
#  define USE_LP_POSITIVITY_LIMITER 0
#endif // USE_LP_POSITIVITY_LIMITER
template <class GV, size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   HatFunctionMomentBasis<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = HatFunctionMomentBasis<double, 1, double, dimRange, 1, 1>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t quad_order = 15;
  static constexpr size_t num_quad_refinements = 0;

#if HAVE_CLP && USE_LP_POSITIVITY_LIMITER
  template <class EigenVectorWrapperType>
  static std::unique_ptr<LpPositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = LpPositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, epsilon);
  }
#else // HAVE_CLP
  template <class EigenVectorWrapperType>
  static std::unique_ptr<PositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = PositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, epsilon);
  }
#endif // HAVE_CLP
};

template <class GV, size_t dimRange, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   PartialMomentBasis<double, 1, double, dimRange, 1, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = PartialMomentBasis<double, 1, double, dimRange, 1, 1>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t quad_order = 15;
  static constexpr size_t num_quad_refinements = 0;

  template <class EigenVectorWrapperType>
  static std::unique_ptr<Dg1dRealizabilityLimitedSlope<GV, double, dimRange, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& basis_functions, const double epsilon)
  {
    using SlopeType = Dg1dRealizabilityLimitedSlope<GV, double, dimRange, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, basis_functions, epsilon);
  }
};

#if HAVE_CLP
template <class GV, size_t order, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   RealSphericalHarmonicsMomentBasis<double, double, order, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = RealSphericalHarmonicsMomentBasis<double, double, order, 3>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t quad_order = 2 * order + 6;
  static constexpr size_t num_quad_refinements = 0;

  template <class EigenVectorWrapperType>
  static std::unique_ptr<LpConvexhullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& basis_functions, const double epsilon)
  {
    using SlopeType = LpConvexhullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, basis_functions, epsilon);
  }
};
#endif

template <class GV, size_t refinements, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   HatFunctionMomentBasis<double, 3, double, refinements, 1, 3>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = HatFunctionMomentBasis<double, 3, double, refinements, 1, 3>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  static constexpr size_t quad_order = 7; // fekete rule number 7
  static constexpr size_t num_quad_refinements = 0;

#if HAVE_CLP && USE_LP_POSITIVITY_LIMITER
  template <class EigenVectorWrapperType>
  static std::unique_ptr<LpPositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = LpPositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, epsilon);
  }
#else // HAVE_CLP
  template <class EigenVectorWrapperType>
  static std::unique_ptr<PositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& /*basis_functions*/, const double epsilon)
  {
    using SlopeType = PositivityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, epsilon);
  }
#endif // HAVE_CLP
};

#if HAVE_QHULL
template <class GV, size_t refinements, class AnalyticalFluxType, class DiscreteFunctionType>
struct RealizabilityLimiterChooser<GV,
                                   PartialMomentBasis<double, 3, double, refinements, 1, 3, 1>,
                                   AnalyticalFluxType,
                                   DiscreteFunctionType>
{
  using MomentBasis = PartialMomentBasis<double, 3, double, refinements, 1, 3>;
  using EntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
  static constexpr size_t quad_order = 3; // fekete rule number 3
  static constexpr size_t num_quad_refinements = 0;

  template <class EigenVectorWrapperType>
  static std::unique_ptr<DgConvexHullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>>
  make_slope(const EntropyFluxType& entropy_flux, const MomentBasis& basis_functions, const double epsilon)
  {
    using SlopeType = DgConvexHullRealizabilityLimitedSlope<GV, MomentBasis, EigenVectorWrapperType>;
    return std::make_unique<SlopeType>(entropy_flux, basis_functions, epsilon);
  }
};
#endif // HAVE_QHULL

// SourceBeam Pn
template <class MomentBasisImp, bool reconstruct>
struct SourceBeamPnExpectedResults;

template <bool reconstruct>
struct SourceBeamPnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33066818456325309 : 0.33107004463414219;
  static constexpr double l2norm = reconstruct ? 0.46157514055648202 : 0.44609169128864046;
  static constexpr double linfnorm = reconstruct ? 1.1553979882432905 : 1.0882801946666183;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33146057542497681 : 0.33146794280840425;
  static constexpr double l2norm = reconstruct ? 0.46411980559363358 : 0.44913032300780292;
  static constexpr double linfnorm = reconstruct ? 0.98904667015384473 : 0.98709215129457029;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct SourceBeamPnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337610927 : 0.33140398337603194;
  static constexpr double l2norm = reconstruct ? 0.47294828933204164 : 0.45667075585121392;
  static constexpr double linfnorm = reconstruct ? 1.0490804598503625 : 0.99004736850989217;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct SourceBeamPnTestCase
{
  using MomentBasis = MomentBasisImp;
  static constexpr size_t dimDomain = MomentBasis::dimDomain;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  using DomainFieldType = typename MomentBasis::DomainFieldType;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using GridType = GridImp;
  using GridViewType = typename GridType::LeafGridView;
  using E = XT::Grid::extract_entity_t<GridViewType>;
  using SpaceType = FiniteVolumeSpace<GridViewType, dimRange, 1, RangeFieldType>;
  using AdvectionSourceSpaceType =
      std::conditional_t<reconstruct, DiscontinuousLagrangeSpace<GridViewType, dimRange, RangeFieldType>, SpaceType>;
  using VectorType = typename Dune::XT::LA::Container<RangeFieldType, Dune::XT::LA::default_backend>::VectorType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ProblemType = SourceBeamPn<E, MomentBasis>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = SourceBeamPnExpectedResults<MomentBasisImp, reconstruction>;
};


// SourceBeam Mn
template <class MomentBasisImp, bool reconstruct>
struct SourceBeamMnExpectedResults;

template <bool reconstruct>
struct SourceBeamMnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140386483109008 : 0.33140386481818457;
  static constexpr double l2norm = reconstruct ? 0.45584415358140384 : 0.44485813650515416;
  static constexpr double linfnorm = reconstruct ? 0.99172157113852111 : 0.98930892893952782;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398330545301 : 0.33140398330533227;
  static constexpr double l2norm = reconstruct ? 0.45584140597017353 : 0.44485191601010715;
  static constexpr double linfnorm = reconstruct ? 0.99172197084890834 : 0.98930925210045084;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct SourceBeamMnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0.33140398337368543 : 0.3314039833756291;
  static constexpr double l2norm = reconstruct ? 0.45583354074069732 : 0.44484887610818585;
  static constexpr double linfnorm = reconstruct ? 0.99172184304625632 : 0.98930905293056492;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct SourceBeamMnTestCase : public SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::GridViewType;
  using ProblemType = SourceBeamMn<GridViewType, MomentBasisImp>;
  using ExpectedResultsType = SourceBeamMnExpectedResults<MomentBasisImp, reconstruct>;
  using RealizabilityLimiterChooserType =
      RealizabilityLimiterChooser<GridViewType, MomentBasisImp, typename ProblemType::FluxType, DiscreteFunctionType>;
};

// PlaneSource Pn
template <class MomentBasisImp, bool reconstruct>
struct PlaneSourcePnExpectedResults;

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.9616518419466558 : 2.7792352623482848;
  static constexpr double linfnorm = reconstruct ? 7.5355813391308644 : 5.9472849007944166;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000240000149;
  static constexpr double l2norm = reconstruct ? 2.8915349919892397 : 2.7676677008555917;
  static constexpr double linfnorm = reconstruct ? 6.9950740716997668 : 5.8904604670932663;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PlaneSourcePnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000239999896 : 2.0000000239999918;
  static constexpr double l2norm = reconstruct ? 2.8799152602279068 : 2.771228836660768;
  static constexpr double linfnorm = reconstruct ? 6.9320887958307775 : 6.0090382693364512;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct PlaneSourcePnTestCase : SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using typename BaseType::E;
  using ProblemType = PlaneSourcePn<E, MomentBasisImp>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PlaneSourcePnExpectedResults<MomentBasisImp, reconstruction>;
};


// PlaneSource Mn
template <class MomentBasisImp, bool reconstruct>
struct PlaneSourceMnExpectedResults;

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<LegendreMomentBasis<double, double, 7>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000240000007 : 2.0000000240000029;
  static constexpr double l2norm = reconstruct ? 2.793246778535933 : 2.746101358507282;
  static constexpr double linfnorm = reconstruct ? 4.9009379502773465 : 5.327698357914608;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<HatFunctionMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = 2.0000000239315696;
  static constexpr double l2norm = reconstruct ? 2.7966600752714887 : 2.7457411547488615;
  static constexpr double linfnorm = reconstruct ? 5.2425259627991894 : 4.9923971272638816;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PlaneSourceMnExpectedResults<PartialMomentBasis<double, 1, double, 8, 1, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 2.0000000239999913 : 2.0000000239999904;
  static constexpr double l2norm = reconstruct ? 2.8215879031834015 : 2.7633864171098814;
  static constexpr double linfnorm = reconstruct ? 6.0674052799351612 : 6.2607864745531092;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct PlaneSourceMnTestCase : SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using typename BaseType::DiscreteFunctionType;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using typename BaseType::GridViewType;
  using ProblemType = PlaneSourceMn<GridViewType, MomentBasisImp>;
  static constexpr RangeFieldType t_end = 0.25;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PlaneSourceMnExpectedResults<MomentBasisImp, reconstruction>;
  using RealizabilityLimiterChooserType =
      RealizabilityLimiterChooser<GridViewType, MomentBasisImp, typename ProblemType::FluxType, DiscreteFunctionType>;
};


// PointSourcePn
template <class MomentBasisImp, bool reconstruct>
struct PointSourcePnExpectedResults
{
  static constexpr double l1norm = 0.;
  static constexpr double l2norm = 0.;
  static constexpr double linfnorm = 0.;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0002337655521842 : 1.0002337655520239;
  static constexpr double l2norm = reconstruct ? 2.6915260598031385 : 2.6808230082176974;
  static constexpr double linfnorm = reconstruct ? 10.353423139916222 : 10.360026140209081;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0002469302942893 : 1.0002469302942893;
  static constexpr double l2norm = reconstruct ? 2.7010206252657687 : 2.6896262652851561;
  static constexpr double linfnorm = reconstruct ? 10.396434946616649 : 10.398129260870316;
  // The matrices in this test case all have eigenvalues [+-0.808311035811965, 0, 0, 0, 0].
  // Thus, the eigenvectors are not unique, and the eigensolvers are extremely sensitive
  // to numerical errors. A difference of 1e-16 in the jacobians entries suffices to
  // result in completely different eigenvectors. In all cases, the eigenvectors are
  // valid eigenvectors to the correct eigenvalues. However, this difference in the
  // eigendecomposition leads to differences in the results with linear reconstruction
  // that are larger than would be expected by pure numerical errors.
  static constexpr double tol = reconstruct ? 1e-5 : 1e-9;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<HatFunctionMomentBasis<double, 3, double, 1, 1, 3>, reconstruct>
{
  // Results with reconstruction not available yet
  static constexpr double l1norm = 1.0002341652578608;
  static constexpr double l2norm = 2.6887381594716606;
  static constexpr double linfnorm = 10.395935217684613;
  // see above
  static constexpr double tol = reconstruct ? 1e-5 : 1e-9;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PartialMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0002469302944583 : 1.0002469302942962;
  static constexpr double l2norm = reconstruct ? 2.6992331534557046 : 2.6888397959991308;
  static constexpr double linfnorm = reconstruct ? 10.39376582567558 : 10.396602604603611;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PointSourcePnExpectedResults<PartialMomentBasis<double, 3, double, 1, 1, 3>, reconstruct>
{
  static_assert(!reconstruct, "Results with reconstruction not available yet!");
  static constexpr double l1norm = 1.0002341652578572;
  static constexpr double l2norm = 2.6888186684779529;
  static constexpr double linfnorm = 10.396523497423775;
  static constexpr double tol = 1e-9;
};


template <class GridImp, class MomentBasisImp, bool reconstruct>
struct PointSourcePnTestCase : SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using typename BaseType::E;
  using ProblemType = PointSourcePn<E, MomentBasisImp>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourcePnExpectedResults<MomentBasisImp, reconstruction>;
};

// CheckerboardPn
template <class MomentBasisImp, bool reconstruct>
struct CheckerboardPnExpectedResults
{
  static constexpr double l1norm = 0.;
  static constexpr double l2norm = 0.;
  static constexpr double linfnorm = 0;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct CheckerboardPnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = 0.35404421488259252;
  static constexpr double l2norm = 0.32914456460176678;
  static constexpr double linfnorm = 0.32888162597935433;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct CheckerboardPnTestCase : SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using typename BaseType::E;
  using ProblemType = CheckerboardPn<E, MomentBasisImp>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = CheckerboardPnExpectedResults<MomentBasisImp, reconstruction>;
};

// ShadowPn
template <class MomentBasisImp, bool reconstruct>
struct ShadowPnExpectedResults
{
  static constexpr double l1norm = 0.;
  static constexpr double l2norm = 0.;
  static constexpr double linfnorm = 0;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct ShadowPnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = 0.5925663469262351;
  static constexpr double l2norm = 0.097668751162279022;
  static constexpr double linfnorm = 0.016484430782606363;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct ShadowPnTestCase : SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamPnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using typename BaseType::E;
  using ProblemType = ShadowPn<E, MomentBasisImp>;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = ShadowPnExpectedResults<MomentBasisImp, reconstruction>;
};


// PointSourceMn
template <class MomentBasisImp, bool reconstruct>
struct PointSourceMnExpectedResults;

template <bool reconstruct>
struct PointSourceMnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0000016440585395 : 1.0000016440583621;
  static constexpr double l2norm = reconstruct ? 2.6896356393271734 : 2.6837985426830504;
  static constexpr double linfnorm = reconstruct ? 10.368884563879304 : 10.375079241728713;
  static constexpr double tol = 1e-9;
};

template <bool reconstruct>
struct PointSourceMnExpectedResults<HatFunctionMomentBasis<double, 3, double, 0, 1, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0000000829624791 : 1.0000000829622864;
  static constexpr double l2norm = reconstruct ? 2.694751941188763 : 2.6892684619955305;
  static constexpr double linfnorm = reconstruct ? 10.379060444346454 : 10.395305896397684;
  // The matrices in this test case all have eigenvalues [+-0.808311035811965, 0, 0, 0, 0].
  // Thus, the eigenvectors are not unique, and the eigensolvers are extremely sensitive
  // to numerical errors. A difference of 1e-16 in the jacobians entries suffices to
  // result in completely different eigenvectors. In all cases, the eigenvectors are
  // valid eigenvectors to the correct eigenvalues. However, this difference in the
  // eigendecomposition leads to differences in the results with linear reconstruction
  // that are larger than would be expected by pure numerical errors.
  static constexpr double tol = reconstruct ? 1e-5 : 1e-9;
};

template <bool reconstruct>
struct PointSourceMnExpectedResults<PartialMomentBasis<double, 3, double, 0, 1, 3, 1>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 1.0000127558028837 : 1.0000127558027134;
  static constexpr double l2norm = reconstruct ? 2.6983746461758269 : 2.6882147627660715;
  static constexpr double linfnorm = reconstruct ? 10.391230462003152 : 10.394186618150204;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct PointSourceMnTestCase : SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using typename BaseType::GridViewType;
  using ProblemType = PointSourceMn<GridViewType, MomentBasisImp>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = PointSourceMnExpectedResults<MomentBasisImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<GridViewType,
                                                                      MomentBasisImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};


// CheckerboardMn
template <class MomentBasisImp, bool reconstruct>
struct CheckerboardMnExpectedResults;

template <bool reconstruct>
struct CheckerboardMnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0. : 0.35404367426148931;
  static constexpr double l2norm = reconstruct ? 0. : 0.32921585327065972;
  static constexpr double linfnorm = reconstruct ? 0. : 0.32895499053937993;
  static constexpr double tol = 1e-9;
};

template <class GridImp, class MomentBasisImp, bool reconstruct>
struct CheckerboardMnTestCase : SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using typename BaseType::GridViewType;
  using ProblemType = CheckerboardMn<GridViewType, MomentBasisImp>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = CheckerboardMnExpectedResults<MomentBasisImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<GridViewType,
                                                                      MomentBasisImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};


// ShadowMn
template <class MomentBasisImp, bool reconstruct>
struct ShadowMnExpectedResults;

template <bool reconstruct>
struct ShadowMnExpectedResults<RealSphericalHarmonicsMomentBasis<double, double, 2, 3>, reconstruct>
{
  static constexpr double l1norm = reconstruct ? 0. : 0.59247415692031025;
  static constexpr double l2norm = reconstruct ? 0. : 0.097642936357459895;
  static constexpr double linfnorm = reconstruct ? 0. : 0.016480937983211638;
  static constexpr double tol = 1e-9;
};


template <class GridImp, class MomentBasisImp, bool reconstruct>
struct ShadowMnTestCase : SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>
{
  using BaseType = SourceBeamMnTestCase<GridImp, MomentBasisImp, reconstruct>;
  using typename BaseType::GridViewType;
  using ProblemType = ShadowMn<GridViewType, MomentBasisImp>;
  using typename BaseType::RangeFieldType;
  static constexpr RangeFieldType t_end = 0.1;
  static constexpr bool reconstruction = reconstruct;
  using ExpectedResultsType = ShadowMnExpectedResults<MomentBasisImp, reconstruction>;
  using RealizabilityLimiterChooserType = RealizabilityLimiterChooser<GridViewType,
                                                                      MomentBasisImp,
                                                                      typename ProblemType::FluxType,
                                                                      typename BaseType::DiscreteFunctionType>;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICTRANSPORT_TESTCASES_HH
