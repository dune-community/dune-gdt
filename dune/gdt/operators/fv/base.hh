// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

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

#include "datahandle.hh"
#include "enums.hh"
#include "realizability.hh"
#include "reconstructed_function.hh"
//#include "reconstruction_specialized.hh"
//#include "reconstruction.hh"
#include "reconstruction_characteristic.hh"
//#include "reconstruction_dense.hh"
#include "slopelimiters.hh"
#include "alphacalculator.hh"

namespace Dune {
namespace GDT {
namespace internal {


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
          size_t reconstruction_order,
          SlopeLimiters slope_lim,
          class RealizabilityLimiterImp>
class AdvectionTraitsBase
{
public:
  static const size_t polOrder = reconstruction_order;
  static const SlopeLimiters slope_limiter = slope_lim;
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef BoundaryValueImp BoundaryValueType;
  typedef BoundaryInfoImp BoundaryInfoType;
  typedef typename std::remove_pointer_t<std::tuple_element_t<0, BoundaryValueType>> DirichletBoundaryValueType;
  typedef RealizabilityLimiterImp RealizabilityLimiterType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename DirichletBoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename DirichletBoundaryValueType::RangeFieldType RangeFieldType;
  typedef RangeFieldType FieldType;
  typedef typename DirichletBoundaryValueType::DomainType DomainType;
  typedef typename AnalyticalFluxType::PartialURangeType JacobianType;
}; // class AdvectionTraitsBase


} // namespace internal


template <class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueImp,
          class BoundaryInfoImp,
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
  //  static_assert(std::is_base_of< ???, BoundaryValueImp >::value,
  //                "BoundaryValueImp has to be derived from ???!");
  //  static_assert(is_discrete_function<SourceImp>::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function<RangeImp>::value, "RangeImp has to be derived from DiscreteFunction!");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef NumericalCouplingFluxImp NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp NumericalBoundaryFluxType;
  typedef BoundaryValueImp BoundaryValueType;
  typedef BoundaryInfoImp BoundaryInfoType;
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
                              const BoundaryInfoType& boundary_info,
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
                               boundary_info,
                               param,
                               std::forward<LocalOperatorArgTypes>(local_operator_args)...)
  {
    this->append(
        local_operator_,
        new XT::Grid::ApplyOn::PartitionSetInnerIntersectionsPrimally<GridLayerType, Dune::Partitions::Interior>());
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
          class RealizabilityLimiterType>
struct AdvectionOperatorApplier
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class BoundaryInfoType,
            class SourceType,
            class RangeType,
            class DomainFieldType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const BoundaryInfoType& boundary_info,
        SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        bool is_linear,
        const Dune::QuadratureRule<DomainFieldType, 1> quadrature1d,
        const Dune::QuadratureRule<DomainFieldType,
                                   std::remove_pointer_t<std::tuple_element_t<0, BoundaryValueType>>::dimDomain - 1>
            intersection_quadrature,
        const bool regularize,
        const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter,
        LocalOperatorArgTypes&&... local_operator_args)
  {
    typedef typename SourceType::SpaceType::GridLayerType GridLayerType;
    typedef typename std::remove_pointer_t<std::tuple_element_t<0, BoundaryValueType>> DirichletBoundaryValueType;
    typedef typename DirichletBoundaryValueType::DomainType DomainType;
    const GridLayerType& grid_layer = source.space().grid_layer();
    const auto& index_set = grid_layer.indexSet();

    if (regularize) {
      auto moment_regularizer = MomentRegularizer<SourceType, AnalyticalFluxType>(source, analytical_flux, param);
      auto walker = XT::Grid::Walker<GridLayerType>(grid_layer);
      walker.append(moment_regularizer);
      walker.walk(true);
      analytical_flux.reset_regularization();
    }

    // evaluate cell averages
    std::vector<typename DirichletBoundaryValueType::RangeType> source_values(index_set.size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto& entity_index = index_set.index(entity);
      const auto& local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    typedef std::
        vector<std::map<DomainType, typename DirichletBoundaryValueType::RangeType, XT::Common::FieldVectorLess>>
            ReconstructedValuesType;
    ReconstructedValuesType reconstructed_values(grid_layer.size(0));

    auto local_reconstruction_operator =
        LocalReconstructionFvOperator<SourceType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>(
            source,
            source_values,
            analytical_flux,
            boundary_values,
            boundary_info,
            param,
            is_linear,
            quadrature1d,
            reconstructed_values);
    auto walker = XT::Grid::Walker<GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);

    // realizability limiting
    if (realizability_limiter) {
      //      realizability_limiter->set_time(param.get("t")[0]);
      //      realizability_limiter->set_regularization_parameters(analytical_flux_.regularization_parameters(),
      //      analytical_flux_.regularization_times());
      realizability_limiter->set_index_set(&index_set);
      realizability_limiter->set_source_values(&source_values);
      realizability_limiter->set_reconstructed_values(&reconstructed_values);
      walker.clear();
      walker.append(*realizability_limiter);
      walker.walk(true);
    }

    if (regularize) {
      auto optimizer = EntropyProblemSolver<AnalyticalFluxType,
                                            SourceType,
                                            DirichletBoundaryValueType::dimDomain,
                                            DirichletBoundaryValueType::dimRange>(
          analytical_flux, source_values, index_set, reconstructed_values, param);
      walker.clear();
      walker.append(optimizer);
      walker.walk(true);
    }

    typedef ReconstructedLocalizableFunction<GridLayerType,
                                             DomainFieldType,
                                             DirichletBoundaryValueType::dimDomain,
                                             typename AnalyticalFluxType::RangeFieldType,
                                             DirichletBoundaryValueType::dimRange>
        ReconstructedLocalizableFunctionType;
    const ReconstructedLocalizableFunctionType reconstructed_function(grid_layer, reconstructed_values);

    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                BoundaryInfoType,
                                ReconstructedLocalizableFunctionType,
                                RangeType>
    localizable_operator(analytical_flux,
                         boundary_values,
                         boundary_info,
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
          class RealizabilityLimiterType>
struct AdvectionOperatorApplier<NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                0,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  template <class AnalyticalFluxType,
            class BoundaryValueType,
            class BoundaryInfoType,
            class SourceType,
            class RangeType,
            class DomainFieldType,
            class... LocalOperatorArgTypes>
  static void
  apply(const AnalyticalFluxType& analytical_flux,
        const BoundaryValueType& boundary_values,
        const BoundaryInfoType& boundary_info,
        const SourceType& source,
        RangeType& range,
        const XT::Common::Parameter& param,
        const bool /*is_linear*/,
        const Dune::QuadratureRule<DomainFieldType, 1> /*quadrature_1d*/,
        const Dune::QuadratureRule<DomainFieldType,
                                   std::remove_pointer_t<std::tuple_element_t<0, BoundaryValueType>>::dimDomain - 1>
            intersection_quadrature,
        const bool /*regularize*/,
        const std::shared_ptr<RealizabilityLimiterType> /*realizability_limiter*/,
        LocalOperatorArgTypes&&... local_operator_args)
  {
    AdvectionLocalizableDefault<AnalyticalFluxType,
                                NumericalCouplingFluxType,
                                NumericalBoundaryFluxType,
                                BoundaryValueType,
                                BoundaryInfoType,
                                SourceType,
                                RangeType>
    localizable_operator(analytical_flux,
                         boundary_values,
                         boundary_info,
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
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  static const size_t dimRangeCols = Traits::dimRangeCols;
  static const size_t polOrder = Traits::polOrder;
  static const SlopeLimiters slope_limiter = Traits::slope_limiter;
  typedef typename Traits::NumericalCouplingFluxType NumericalCouplingFluxType;
  typedef typename Traits::NumericalBoundaryFluxType NumericalBoundaryFluxType;
  typedef typename Traits::RealizabilityLimiterType RealizabilityLimiterType;

  typedef Dune::QuadratureRule<DomainFieldType, 1> OnedQuadratureType;
  typedef Dune::QuadratureRule<DomainFieldType, dimDomain - 1> IntersectionQuadratureType;

public:
  AdvectionOperatorBase(const AnalyticalFluxType& analytical_flux,
                        const BoundaryValueType& boundary_values,
                        const BoundaryInfoType& boundary_info,
                        const bool is_linear)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , boundary_info_(boundary_info)
    , is_linear_(is_linear)
    , quadrature_1d_(default_quadrature_helper<>::get())
    , intersection_quadrature_(quadrature_helper<>::get(quadrature_1d_))
    , regularize_(false)
  {
  }

  AdvectionOperatorBase(const AnalyticalFluxType& analytical_flux,
                        const BoundaryValueType& boundary_values,
                        const BoundaryInfoType& boundary_info,
                        const bool is_linear,
                        const OnedQuadratureType& quadrature_1d,
                        const bool regularize,
                        const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , boundary_info_(boundary_info)
    , is_linear_(is_linear)
    , quadrature_1d_(quadrature_1d)
    , intersection_quadrature_(quadrature_helper<>::get(quadrature_1d_))
    , regularize_(regularize)
    , realizability_limiter_(realizability_limiter)
  {
  }

  template <class SourceType, class RangeType, class... Args>
  void apply(SourceType& source, RangeType& range, const XT::Common::Parameter& param, Args&&... args) const
  {
    internal::AdvectionOperatorApplier<NumericalCouplingFluxType,
                                       NumericalBoundaryFluxType,
                                       polOrder,
                                       slope_limiter,
                                       RealizabilityLimiterType>::apply(analytical_flux_,
                                                                        boundary_values_,
                                                                        boundary_info_,
                                                                        source,
                                                                        range,
                                                                        param,
                                                                        is_linear_,
                                                                        quadrature_1d_,
                                                                        intersection_quadrature_,
                                                                        regularize_,
                                                                        realizability_limiter_,
                                                                        std::forward<Args>(args)...);
  }

  void set_1d_quadrature(const OnedQuadratureType& quadrature)
  {
    quadrature_1d_ = quadrature;
    intersection_quadrature_ = quadrature_helper<>::get(quadrature);
  }

  void set_realizability_limiter(const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter)
  {
    realizability_limiter_ = realizability_limiter;
  }

  static OnedQuadratureType default_1d_quadrature()
  {
    return default_quadrature_helper<>::get();
  }

  const AnalyticalFluxType& analytical_flux() const
  {
    return analytical_flux_;
  }

private:
  template <size_t reconstructionOrder = polOrder, class anything = void>
  struct default_quadrature_helper
  {
    static OnedQuadratureType get()
    {
      return Dune::QuadratureRules<DomainFieldType, 1>::rule(Dune::GeometryType(Dune::GeometryType::BasicType::cube, 1),
                                                             2 * polOrder);
    }
  };

  template <class anything>
  struct default_quadrature_helper<1, anything>
  {
    static OnedQuadratureType get()
    {
      OnedQuadratureType quadrature;
      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5, 1.));
      //      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. - 1. / std::sqrt(3)), 0.5));
      //      quadrature.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0.5 * (1. + 1. / std::sqrt(3)), 0.5));
      return quadrature;
    }
  };

  template <size_t domainDim = dimDomain, class anything = void>
  struct quadrature_helper;

  template <class anything>
  struct quadrature_helper<1, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1> get(const OnedQuadratureType& /*quadrature_1d*/)
    {
      Dune::QuadratureRule<DomainFieldType, dimDomain - 1> ret;
      ret.push_back(Dune::QuadraturePoint<DomainFieldType, 0>(FieldVector<DomainFieldType, 0>(0), 1));
      return ret;
    }
  };

  template <class anything>
  struct quadrature_helper<2, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1> get(const OnedQuadratureType& quadrature_1d)
    {
      return quadrature_1d;
    }
  };

  template <class anything>
  struct quadrature_helper<3, anything>
  {
    static Dune::QuadratureRule<DomainFieldType, dimDomain - 1> get(const OnedQuadratureType& quadrature_1d)
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
  const BoundaryInfoType& boundary_info_;
  const bool is_linear_;
  OnedQuadratureType quadrature_1d_;
  IntersectionQuadratureType intersection_quadrature_;
  const bool regularize_;
  std::shared_ptr<RealizabilityLimiterType> realizability_limiter_;
}; // class AdvectionOperatorBase<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_BASE_HH
