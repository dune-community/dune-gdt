// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_OPERATORS_FV_BASE_HH
#define DUNE_GDT_OPERATORS_FV_BASE_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/local/operators/fv.hh>
#include <dune/gdt/operators/base.hh>

#include "reconstructed_function.hh"
#include "reconstruction.hh"
#include "slopelimiters.hh"
#include "realizability.hh"

namespace Dune {
namespace GDT {


enum class NumericalFluxes
{
  godunov,
  godunov_with_reconstruction,
  laxfriedrichs,
  laxfriedrichs_with_reconstruction,
  local_laxfriedrichs,
  local_laxfriedrichs_with_reconstruction,
  kinetic
};


namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          size_t reconstructionOrder,
          SlopeLimiters slope_lim,
          bool realizability_lim,
          class BasisFunctionImp,
          class EigenSolverImp>
class AdvectionTraitsBase
{
public:
  static const size_t polOrder = reconstructionOrder;
  static const SlopeLimiters slope_limiter = slope_lim;
  static const bool realizability_limiting = realizability_lim;
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueType;
  typedef BasisFunctionImp BasisFunctionType;
  typedef EigenSolverImp EigenSolverType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename BoundaryValueFunctionImp::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueFunctionImp::RangeFieldType RangeFieldType;
  typedef RangeFieldType FieldType;
  typedef typename BoundaryValueFunctionImp::DomainType DomainType;
  typedef typename AnalyticalFluxType::PartialURangeType JacobianType;
}; // class AdvectionTraitsBase


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueFunctionImp,
          class SourceImp,
          class RangeImp>
class AdvectionLocalizableDefault
    : public Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp>
{
  typedef Dune::GDT::LocalizableOperatorBase<typename RangeImp::SpaceType::GridLayerType, SourceImp, RangeImp> BaseType;

  static_assert(is_local_numerical_coupling_flux<NumericalCouplingFluxImp>::value,
                "NumericalCouplingFluxImp has to be derived from LocalNumericalCouplingFluxInterface!");
  static_assert(is_local_numerical_boundary_flux<NumericalBoundaryFluxImp>::value,
                "NumericalBoundaryFluxImp has to be derived from LocalNumericalBoundaryFluxInterface!");
  //  static_assert(std::is_base_of< ???, BoundaryValueFunctionImp >::value,
  //                "BoundaryValueFunctionImp has to be derived from ???!");
  //  static_assert(is_discrete_function<SourceImp>::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
  typedef typename SourceType::RangeFieldType RangeFieldType;
  typedef typename RangeType::SpaceType::GridLayerType GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename Dune::GDT::LocalCouplingFvOperator<NumericalCouplingFluxType> LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFvOperator<NumericalBoundaryFluxType> LocalBoundaryOperatorType;

  template <class QuadratureRuleType, class... LocalOperatorArgTypes>
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const BoundaryValueType& boundary_values,
                              const SourceType& source,
                              RangeType& range,
                              const XT::Common::Parameter& param,
                              const QuadratureRuleType& quadrature_rule,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(range.space().grid_layer(), source, range)
    , local_operator_(
          quadrature_rule, analytical_flux, param, std::forward<LocalOperatorArgTypes>(local_operator_args)...)
    , local_boundary_operator_(quadrature_rule,
                               analytical_flux,
                               boundary_values,
                               param,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...)
  {
    this->append(local_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridLayerType>());
    this->append(local_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GridLayerType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<GridLayerType>());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


namespace internal {


template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          size_t polOrder,
          SlopeLimiters slope_limiter,
          bool realizability_limiting,
          class EigenSolverType>
struct AdvectionOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class SourceType,
            class RangeType,
            class RangeFieldType,
            class DomainFieldType,
            class BasisFunctionType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        bool is_linear,
        const Dune::QuadratureRule<DomainFieldType, 1> intersection_quadrature_1d,
        const Dune::QuadratureRule<DomainFieldType, BoundaryValueType::dimDomain - 1> intersection_quadrature,
        const Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain>& quadrature,
        const RangeFieldType epsilon,
        const std::shared_ptr<const BasisFunctionType> basis_functions,
        LocalOperatorArgTypes&&... local_operator_args)

  {
    typedef typename SourceType::SpaceType::GridLayerType GridLayerType;
    typedef typename BoundaryValueType::DomainType DomainType;
    static const size_t dimDomain = BoundaryValueType::dimDomain;
    static const size_t dimRange = BoundaryValueType::dimRange;
    const GridLayerType& grid_layer = source.space().grid_layer();

    // evaluate cell averages
    std::vector<typename BoundaryValueType::RangeType> source_values(grid_layer.indexSet().size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto& entity_index = grid_layer.indexSet().index(entity);
      const auto& local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    std::vector<std::map<DomainType, typename BoundaryValueType::RangeType, XT::Common::FieldVectorLess>>
        reconstructed_values(grid_layer.size(0));

    auto local_reconstruction_operator = LocalReconstructionFvOperator<GridLayerType,
                                                                       AnalyticalFluxType,
                                                                       BoundaryValueType,
                                                                       polOrder,
                                                                       slope_limiter,
                                                                       EigenSolverType>(source_values,
                                                                                        analytical_flux,
                                                                                        boundary_values,
                                                                                        grid_layer,
                                                                                        param,
                                                                                        is_linear,
                                                                                        intersection_quadrature_1d,
                                                                                        reconstructed_values);
    auto walker = XT::Grid::Walker<GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);

    if (realizability_limiting) {
      assert(basis_functions);
      // do limiting for realizability in M_N models
      //      auto local_realizability_limiter =
      //          LocalRealizabilityLimiter<SourceType, BasisFunctionType, BasisFunctionType::dimDomain, dimRange>(
      //              source, reconstructed_values, *basis_functions, quadrature, epsilon);
      auto local_realizability_limiter =
          LocalRealizabilityLimiterLP<SourceType, BasisFunctionType, BasisFunctionType::dimDomain, dimRange>(
              source, reconstructed_values, *basis_functions, quadrature, epsilon);
      walker.clear();
      walker.append(local_realizability_limiter);
      walker.walk(true);
    }

    typedef ReconstructedLocalizableFunction<GridLayerType, RangeFieldType, dimDomain, RangeFieldType, dimRange>
        ReconstructedLocalizableFunctionType;
    const ReconstructedLocalizableFunctionType reconstructed_function(grid_layer, reconstructed_values);

    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                ReconstructedLocalizableFunctionType,
                                RangeType>
        localizable_operator(analytical_flux,
                             boundary_values,
                             reconstructed_function,
                             range,
                             param,
                             intersection_quadrature,
                             std::forward<LocalOperatorArgTypes>(local_operator_args)...);
    localizable_operator.apply(true);
  }
}; // struct AdvectionOperatorApplier<..., polOrder>0,...>

template <class NumericalCouplingFluxType,
          class NumericalBoundaryFluxType,
          SlopeLimiters slope_limiter,
          bool realizability_limiting,
          class EigenSolverType>
struct AdvectionOperatorApplier<NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                0,
                                slope_limiter,
                                realizability_limiting,
                                EigenSolverType>
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class SourceType,
            class RangeType,
            class RangeFieldType,
            class DomainFieldType,
            class BasisFunctionType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        const bool /*is_linear*/,
        const Dune::QuadratureRule<DomainFieldType, 1> /*intersection_quadrature_1d*/,
        const Dune::QuadratureRule<DomainFieldType, BoundaryValueType::dimDomain - 1> intersection_quadrature,
        const Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain> /*quadrature*/,
        const RangeFieldType /*epsilon*/,
        const std::shared_ptr<const BasisFunctionType> /*basis_functions*/,
        LocalOperatorArgTypes&&... local_operator_args)
  {
    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                SourceType,
                                RangeType>
        localizable_operator(analytical_flux,
                             boundary_values,
                             source,
                             range,
                             param,
                             intersection_quadrature,
                             std::forward<LocalOperatorArgTypes>(local_operator_args)...);
    localizable_operator.apply(true);
  }
}; // struct AdvectionOperatorApplier<..., polOrder=0,...>


} // namespace internal


template <class Traits>
class AdvectionOperatorBase
{
public:
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const size_t polOrder = Traits::polOrder;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  static const bool realizability_limiting = Traits::realizability_limiting;
  typedef typename Traits::BasisFunctionType BasisFunctionType;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;
  typedef typename Traits::EigenSolverType EigenSolverType;

  typedef Dune::QuadratureRule<DomainFieldType, 1> Intersection1dQuadratureType;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain - 1> IntersectionQuadratureType;
  typedef Dune::QuadratureRule<DomainFieldType, BasisFunctionType::dimDomain> QuadratureType;

public:
  AdvectionOperatorBase(const AnalyticalFluxType& analytical_flux,
                        const BoundaryValueType& boundary_values,
                        const bool is_linear)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
    , intersection_1d_quadrature_(helper<>::default_1d_quadrature())
    , intersection_quadrature_(helper2<>::get_quadrature(intersection_1d_quadrature_))
    , epsilon_(1e-14)
  {
  }

  template <class SourceType, class RangeType, class... Args>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param, Args&&... args) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       polOrder,
                                       slope_limiter,
                                       realizability_limiting,
                                       EigenSolverType>::apply(analytical_flux_,
                                                               boundary_values_,
                                                               source,
                                                               range,
                                                               param,
                                                               is_linear_,
                                                               intersection_1d_quadrature_,
                                                               intersection_quadrature_,
                                                               quadrature_,
                                                               epsilon_,
                                                               basis_functions_,
                                                               std::forward<Args>(args)...);
  }

  void set_1d_quadrature(const Intersection1dQuadratureType& quadrature)
  {
    intersection_1d_quadrature_ = quadrature;
    intersection_quadrature_ = helper2<>::get_quadrature(quadrature);
  }

  void set_quadrature(const QuadratureType& quadrature)
  {
    quadrature_ = quadrature;
  }

  void set_epsilon(const RangeFieldType& epsilon)
  {
    epsilon_ = epsilon;
  }

  void set_basisfunctions(const std::shared_ptr<const BasisFunctionType> basis_functions)
  {
    basis_functions_ = basis_functions;
  }

private:
  template <size_t reconstructionOrder = polOrder, class anything = void>
  struct helper
  {
    static Intersection1dQuadratureType default_1d_quadrature()
    {
      return Dune::QuadratureRules<DomainFieldType, 1>::rule(Dune::GeometryType(Dune::GeometryType::BasicType::cube, 1),
                                                             2 * polOrder);
    }
  };

  template <class anything>
  struct helper<1, anything>
  {
    static Intersection1dQuadratureType default_1d_quadrature()
    {
      // get 1D quadrature rules
      Intersection1dQuadratureType quadrature;
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. - 1. / std::sqrt(3)), 0.5));
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. + 1. / std::sqrt(3)), 0.5));
      return quadrature;
    }
  };

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper2;

  template <class anything>
  struct helper2<1, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& /*quadrature_1d*/)
    {
      Dune::QuadratureRule<DomainFieldType, dimDomain - 1> ret;
      ret.push_back(Dune::QuadraturePoint<DomainFieldType, 0>(FieldVector<DomainFieldType, 0>(0), 1));
      return ret;
    }
  };

  template <class anything>
  struct helper2<2, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& quadrature_1d)
    {
      return quadrature_1d;
    }
  };

  template <class anything>
  struct helper2<3, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1>
    get_quadrature(const Intersection1dQuadratureType& quadrature_1d)
    {
      Dune::QuadratureRule<DomainFieldType, dimDomain - 1> ret;
      for (size_t ii = 0; ii < quadrature_1d.size(); ++ii)
        for (size_t jj = 0; jj < quadrature_1d.size(); ++jj)
          ret.push_back(Dune::QuadraturePoint<DomainFieldType, dimDomain - 1>(
              {quadrature_1d[ii].position()[0], quadrature_1d[jj].position()[0]},
              quadrature_1d[ii].weight() * quadrature_1d[jj].weight()));
      return ret;
    }
  };

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const bool is_linear_;
  Intersection1dQuadratureType intersection_1d_quadrature_;
  IntersectionQuadratureType intersection_quadrature_;
  QuadratureType quadrature_;
  std::shared_ptr<const BasisFunctionType> basis_functions_;
  RangeFieldType epsilon_;
}; // class AdvectionOperatorBase<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_BASE_HH
