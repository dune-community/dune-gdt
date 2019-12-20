// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_MOMENTMODELS_DENSITYEVALUATOR_HH
#define DUNE_GDT_MOMENTMODELS_DENSITYEVALUATOR_HH

#include <string>
#include <functional>

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class MomentBasis, SlopeLimiterType slope>
class LocalDensityEvaluator : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using BaseType = XT::Grid::ElementFunctor<typename SpaceType::GridViewType>;

public:
  using GridViewType = typename SpaceType::GridViewType;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GridViewType, MomentBasis, slope>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  static const size_t dimFlux = EntropyFluxType::dimFlux;
  static const size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using DomainType = FieldVector<RangeFieldType, dimFlux>;
  using BoundaryDistributionType = std::function<std::function<RangeFieldType(const DomainType&)>(const DomainType&)>;

  explicit LocalDensityEvaluator(const SpaceType& space,
                                 const VectorType& alpha_dofs,
                                 VectorType& range_dofs,
                                 EntropyFluxType& analytical_flux,
                                 const BoundaryDistributionType& boundary_distribution,
                                 const RangeFieldType min_acceptable_density,
                                 const XT::Common::Parameter& param)
    : space_(space)
    , alpha_(space_, alpha_dofs, "alpha")
    , range_(space_, range_dofs, "regularized alpha")
    , local_alpha_(alpha_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , boundary_distribution_(boundary_distribution)
    , min_acceptable_density_(min_acceptable_density)
    , param_(param)
    , index_set_(space_.grid_view().indexSet())
  {}

  explicit LocalDensityEvaluator(LocalDensityEvaluator& other)
    : BaseType(other)
    , space_(other.space_)
    , alpha_(space_, other.alpha_.dofs().vector(), "source")
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_alpha_(alpha_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , boundary_distribution_(other.boundary_distribution_)
    , min_acceptable_density_(other.min_acceptable_density_)
    , param_(other.param_)
    , index_set_(space_.grid_view().indexSet())
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalDensityEvaluator(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_alpha_->bind(entity);
    local_range_->bind(entity);
    const auto entity_index = index_set_.index(entity);
    const auto& local_alpha_dofs = local_alpha_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      alpha_tmp_[ii] = local_alpha_dofs.get_entry(ii);
    analytical_flux_.store_evaluations(entity_index, alpha_tmp_, min_acceptable_density_);
    auto& local_range_dofs = local_range_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_dofs.set_entry(ii, alpha_tmp_[ii]);
    for (auto&& intersection : Dune::intersections(space_.grid_view(), entity))
      if (intersection.boundary())
        analytical_flux_.store_boundary_evaluations(
            boundary_distribution_(intersection.geometry().center()), entity_index, intersection.indexInInside());
    analytical_flux_.set_eta_ast_pointers();
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType alpha_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_alpha_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  EntropyFluxType& analytical_flux_;
  const BoundaryDistributionType& boundary_distribution_;
  const RangeFieldType min_acceptable_density_;
  const XT::Common::Parameter& param_;
  const typename SpaceType::GridViewType::IndexSet& index_set_;
  XT::Common::FieldVector<RangeFieldType, dimRange> alpha_tmp_;
}; // class LocalDensityEvaluator<...>

template <class MomentBasisImp,
          class SpaceImp,
          SlopeLimiterType slope,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class DensityEvaluator
  : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>
{
  using BaseType = OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>;

public:
  using typename BaseType::VectorType;
  using MomentBasis = MomentBasisImp;
  using SpaceType = SpaceImp;
  using SourceSpaceType = SpaceImp;
  using RangeSpaceType = SpaceImp;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<typename SpaceType::GridViewType, MomentBasis, slope>;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using LocalDensityEvaluatorType = LocalDensityEvaluator<SpaceType, VectorType, MomentBasis, slope>;
  using BoundaryDistributionType = typename LocalDensityEvaluatorType::BoundaryDistributionType;

  DensityEvaluator(EntropyFluxType& analytical_flux,
                   const SpaceType& space,
                   const BoundaryDistributionType& boundary_distribution,
                   const RangeFieldType min_acceptable_density)
    : analytical_flux_(analytical_flux)
    , space_(space)
    , boundary_distribution_(boundary_distribution)
    , min_acceptable_density_(min_acceptable_density)
  {
    analytical_flux_.prepare_storage(space_.grid_view());
  }

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  void apply(const VectorType& alpha, VectorType& range, const XT::Common::Parameter& param = {}) const override final
  {
    LocalDensityEvaluatorType local_density_evaluator(
        space_, alpha, range, analytical_flux_, boundary_distribution_, min_acceptable_density_, param);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_density_evaluator);
    walker.walk(true);
  } // void apply(...)

private:
  EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
  const BoundaryDistributionType& boundary_distribution_;
  const RangeFieldType min_acceptable_density_;
}; // class DensityEvaluator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_DENSITYEVALUATOR_HH
