// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_KINETIC_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_KINETIC_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker.hh>

#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/tools/discretevalued-grid-function.hh>

#include "slopes.hh"
#include "internal.hh"

namespace Dune {
namespace GDT {


template <class GV, class AnalyticalFluxType, class BoundaryValueType, class LocalVectorType>
class LocalPointwiseLinearKineticReconstructionOperator : public XT::Grid::ElementFunctor<GV>
{
  using ThisType = LocalPointwiseLinearKineticReconstructionOperator;
  using BaseType = XT::Grid::ElementFunctor<GV>;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using EntityType = typename GV::template Codim<0>::Entity;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using StencilType = DynamicVector<LocalVectorType>;
  using StencilsType = std::vector<StencilType>;
  using DomainType = typename BoundaryValueType::DomainType;
  using RangeType = typename BoundaryValueType::RangeReturnType;
  static constexpr size_t stencil_size = 3;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, RangeFieldType>;

public:
  explicit LocalPointwiseLinearKineticReconstructionOperator(ReconstructedFunctionType& reconstructed_function,
                                                             const GV& grid_view,
                                                             const AnalyticalFluxType& analytical_flux,
                                                             const std::vector<LocalVectorType>& source_values,
                                                             const BoundaryValueType& boundary_values,
                                                             const XT::Common::Parameter& param)
    : reconstructed_function_(reconstructed_function)
    , grid_view_(grid_view)
    , analytical_flux_(analytical_flux)
    , source_values_(source_values)
    , boundary_values_(boundary_values)
    , param_(param)
    , stencils_(dimDomain, StencilType(stencil_size))
  {}

  LocalPointwiseLinearKineticReconstructionOperator(const ThisType& other)
    : BaseType(other)
    , reconstructed_function_(other.reconstructed_function_)
    , grid_view_(other.grid_view_)
    , analytical_flux_(other.analytical_flux_)
    , source_values_(other.source_values_)
    , boundary_values_(other.boundary_values_)
    , param_(other.param_)
    , stencils_(dimDomain, StencilType(stencil_size))
  {}

  virtual XT::Grid::ElementFunctor<GV>* copy() override final
  {
    return new ThisType(*this);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!fill_stencils(entity))
      return;

    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      analytical_flux_.calculate_minmod_density_reconstruction(stencils_[dd], local_reconstructed_values, dd);
  } // void apply_local(...)

private:
  bool fill_stencils(const EntityType& entity)
  {
    const auto entity_index = grid_view_.indexSet().index(entity);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      stencils_[dd][1] = source_values_[entity_index];
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const size_t dd = intersection.indexInInside() / 2;
      const size_t index = (intersection.indexInInside() % 2) * 2;
      if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
        stencils_[dd][index] = boundary_values_.evaluate(intersection.geometry().center());
      } else if (intersection.neighbor()) { // inner and periodic intersections {
        stencils_[dd][index] = source_values_[grid_view_.indexSet().index(intersection.outside())];
      } else if (!intersection.neighbor() && !intersection.boundary()) { // processor boundary
        return false;
      } else {
        DUNE_THROW(Dune::NotImplemented, "This should not happen!");
      }
    } // intersections
    return true;
  } // void fill_stencils(...)

  ReconstructedFunctionType& reconstructed_function_;
  const GV& grid_view_;
  const AnalyticalFluxType& analytical_flux_;
  const std::vector<LocalVectorType>& source_values_;
  const BoundaryValueType& boundary_values_;
  const XT::Common::Parameter& param_;
  StencilsType stencils_;
}; // class LocalPointwiseLinearKineticReconstructionOperator

// Does not reconstruct a full first-order DG function, but only stores the reconstructed values at the intersection
// centers. This avoids the interpolation in this operator and the evaluation of the reconstructed function in the
// finite volume operator which are both quite expensive for large dimRange.
template <class GV, class BoundaryValueImp, class AnalyticalFluxType, class VectorType, class LocalVectorType>
class PointwiseLinearKineticReconstructionOperator
{
public:
  using BoundaryValueType = BoundaryValueImp;
  using E = XT::Grid::extract_entity_t<GV>;
  using EigenVectorWrapperType = DummyEigenVectorWrapper<LocalVectorType>;
  using R = typename BoundaryValueType::R;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using SpaceType = SpaceInterface<GV, dimRange, 1, R>;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, R>;
  using ReconstructedValuesType = std::vector<typename ReconstructedFunctionType::LocalFunctionValuesType>;

  PointwiseLinearKineticReconstructionOperator(const BoundaryValueType& boundary_values,
                                               const SpaceType& space,
                                               const AnalyticalFluxType& analytical_flux)
    : boundary_values_(boundary_values)
    , space_(space)
    , analytical_flux_(analytical_flux)
  {}

  bool linear() const
  {
    return false;
  }

  const SpaceType& space() const
  {
    return space_;
  }

  void apply(const VectorType& source, ReconstructedFunctionType& range, const XT::Common::Parameter& param) const
  {
    // evaluate cell averages
    const auto& grid_view = space_.grid_view();
    const auto& index_set = grid_view.indexSet();
    std::vector<LocalVectorType> source_values(index_set.size(0));
    // Contains the evaluations of exp(alpha^T b(v_i)) at each quadrature point v_i
    auto source_func = make_discrete_function(space_, source);
    const auto local_source = source_func.local_function();
    for (const auto& entity : Dune::elements(grid_view)) {
      local_source->bind(entity);
      const auto entity_index = index_set.index(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator =
        LocalPointwiseLinearKineticReconstructionOperator<GV, AnalyticalFluxType, BoundaryValueType, LocalVectorType>(
            range, grid_view, analytical_flux_, source_values, boundary_values_, param);
    auto walker = XT::Grid::Walker<GV>(grid_view);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  const BoundaryValueType& boundary_values_;
  const SpaceType& space_;
  const AnalyticalFluxType& analytical_flux_;
}; // class PointwiseLinearKineticReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_KINETIC_HH
