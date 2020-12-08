// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_KINETIC_HH
#define DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_KINETIC_HH

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


template <class GV, class AnalyticalFluxType, class LocalVectorType>
class LocalPointwiseLinearKineticReconstructionOperator : public XT::Grid::ElementFunctor<GV>
{
  using ThisType = LocalPointwiseLinearKineticReconstructionOperator;
  using BaseType = XT::Grid::ElementFunctor<GV>;
  using IndexSetType = typename GV::IndexSet;
  static constexpr size_t dimDomain = GV::dimension;
  static constexpr size_t dimRange = AnalyticalFluxType::state_dim;
  using EntityType = typename GV::template Codim<0>::Entity;
  static constexpr size_t stencil_size = 3;
  using StencilType = XT::Common::FieldVector<size_t, stencil_size>;
  using StencilsType = FieldVector<StencilType, dimDomain>;
  using RangeFieldType = typename AnalyticalFluxType::RangeFieldType;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, RangeFieldType>;

public:
  explicit LocalPointwiseLinearKineticReconstructionOperator(ReconstructedFunctionType& reconstructed_function,
                                                             const GV& grid_view,
                                                             const AnalyticalFluxType& analytical_flux,
                                                             const XT::Common::Parameter& param)
    : reconstructed_function_(reconstructed_function)
    , grid_view_(grid_view)
    , index_set_(grid_view_.indexSet())
    , analytical_flux_(analytical_flux)
    , param_(param)
  {}

  LocalPointwiseLinearKineticReconstructionOperator(const ThisType& other)
    : BaseType(other)
    , reconstructed_function_(other.reconstructed_function_)
    , grid_view_(other.grid_view_)
    , index_set_(grid_view_.indexSet())
    , analytical_flux_(other.analytical_flux_)
    , param_(other.param_)
  {}

  XT::Grid::ElementFunctor<GV>* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!fill_stencils(entity))
      return;

    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      analytical_flux_.calculate_reconstructed_fluxes(
          stencils_[dd], boundary_dirs_[dd], local_reconstructed_values, dd);
  } // void apply_local(...)

private:
  bool fill_stencils(const EntityType& entity)
  {
    for (size_t dd = 0; dd < dimDomain; ++dd)
      stencils_[dd][1] = index_set_.index(entity);
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const size_t dd = intersection.indexInInside() / 2;
      const size_t index = (intersection.indexInInside() % 2) * 2;
      if (intersection.boundary() && !intersection.neighbor()) { // boundary intersections
        stencils_[dd][index] = size_t(-1);
        boundary_dirs_[dd][index] = intersection.indexInInside() % 2;
      } else if (intersection.neighbor()) { // inner and periodic intersections {
        stencils_[dd][index] = index_set_.index(intersection.outside());
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
  const IndexSetType& index_set_;
  const AnalyticalFluxType& analytical_flux_;
  const XT::Common::Parameter& param_;
  StencilsType stencils_;
  XT::Common::FieldVector<XT::Common::FieldVector<size_t, stencil_size>, dimDomain> boundary_dirs_;

}; // class LocalPointwiseLinearKineticReconstructionOperator


// Does not reconstruct a full first-order DG function, but only stores the reconstructed values at the intersection
// centers. This avoids the interpolation in this operator and the evaluation of the reconstructed function in the
// finite volume operator which are both quite expensive for large dimRange.
template <class GV, class AnalyticalFluxType, class VectorType, class LocalVectorType>
class PointwiseLinearKineticReconstructionOperator
{
public:
  using E = XT::Grid::extract_entity_t<GV>;
  using EigenVectorWrapperType = internal::DummyEigenVectorWrapper<AnalyticalFluxType, LocalVectorType>;
  using R = typename AnalyticalFluxType::R;
  static constexpr size_t dimDomain = GV::dimension;
  static constexpr size_t dimRange = AnalyticalFluxType::state_dim;
  using SpaceType = SpaceInterface<GV, dimRange, 1, R>;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, R>;
  using ReconstructedValuesType = std::vector<typename ReconstructedFunctionType::LocalFunctionValuesType>;

  PointwiseLinearKineticReconstructionOperator(const SpaceType& space, const AnalyticalFluxType& analytical_flux)
    : space_(space)
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

  void apply(const VectorType& /*source*/, ReconstructedFunctionType& range, const XT::Common::Parameter& param) const
  {
    // do reconstruction
    const auto& grid_view = space_.grid_view();
    auto local_reconstruction_operator =
        LocalPointwiseLinearKineticReconstructionOperator<GV, AnalyticalFluxType, LocalVectorType>(
            range, grid_view, analytical_flux_, param);
    auto walker = XT::Grid::Walker<GV>(grid_view);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  const SpaceType& space_;
  const AnalyticalFluxType& analytical_flux_;
}; // class PointwiseLinearKineticReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_KINETIC_HH
