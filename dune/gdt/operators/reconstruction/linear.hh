// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_HH

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


template <class AnalyticalFluxType, class BoundaryValueType, class GV, class EigenvectorWrapperType>
class LinearSlopeElementFunctor : public XT::Grid::ElementFunctor<GV>
{
  using ThisType = LinearSlopeElementFunctor;
  using BaseType = XT::Grid::ElementFunctor<GV>;

public:
  using typename BaseType::E;
  using SlopeType = SlopeBase<E, EigenvectorWrapperType>;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using StencilType = DynamicVector<LocalVectorType>;
  using StencilsType = std::vector<StencilType>;
  using DomainType = typename BoundaryValueType::DomainType;
  using RangeType = typename BoundaryValueType::RangeReturnType;
  static constexpr size_t d = BoundaryValueType::d;
  static constexpr size_t stencil_size = 3;

  LinearSlopeElementFunctor(const GV& grid_view,
                            const std::vector<LocalVectorType>& source_values,
                            const BoundaryValueType& boundary_values,
                            const AnalyticalFluxType& analytical_flux,
                            const SlopeType& slope,
                            const XT::Common::Parameter& param,
                            const bool flux_is_affine = false)

    : grid_view_(grid_view)
    , source_values_(source_values)
    , boundary_values_(boundary_values)
    , analytical_flux_(analytical_flux)
    , slope_(slope.copy())
    , param_(param)
    , eigenvector_wrapper_(analytical_flux, flux_is_affine)
    , slopes_(d)
    , stencils_(d, StencilType(stencil_size))
  {}

  LinearSlopeElementFunctor(const LinearSlopeElementFunctor& other)
    : BaseType(other)
    , grid_view_(other.grid_view_)
    , source_values_(other.source_values_)
    , boundary_values_(other.boundary_values_)
    , analytical_flux_(other.analytical_flux_)
    , slope_(other.slope_->copy())
    , param_(other.param_)
    , eigenvector_wrapper_(analytical_flux_, other.eigenvector_wrapper_.affine())
    , slopes_(d)
    , stencils_(d, StencilType(stencil_size))
  {}

  XT::Grid::ElementFunctor<GV>* copy() override final
  {
    return new ThisType(*this);
  }

  ThisType* copy_derived()
  {
    return new ThisType(*this);
  }

  void apply_local(const E& entity) override final
  {
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!fill_stencils(entity))
      return;

    // get eigenvectors of flux
    if (analytical_flux_.x_dependent())
      x_local_ = entity.geometry().local(entity.geometry().center());
    const auto entity_index = grid_view_.indexSet().index(entity);
    eigenvector_wrapper_.compute_eigenvectors(entity, x_local_, source_values_[entity_index], param_);

    for (size_t dd = 0; dd < d; ++dd) {
      // no need to reconstruct in all directions, as we are only regarding the center of the face, which will
      // always have the same value assigned, independent of the slope in the other directions
      slopes_[dd] = slope_->get(entity, stencils_[dd], eigenvector_wrapper_, dd);
    } // dd
  }

  const std::vector<RangeType>& slopes() const
  {
    return slopes_;
  }

  const LocalVectorType& u_entity() const
  {
    return stencils_[0][1];
  }

private:
  bool fill_stencils(const E& entity)
  {
    const auto entity_index = grid_view_.indexSet().index(entity);
    for (size_t dd = 0; dd < d; ++dd)
      stencils_[dd][1] = source_values_[entity_index];
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const size_t dd = intersection.indexInInside() / 2;
      const size_t index = (intersection.indexInInside() % 2) * 2;
      if (intersection.boundary() && !intersection.neighbor()) // boundary intersections
        stencils_[dd][index] = boundary_values_.evaluate(intersection.geometry().center());
      else if (intersection.neighbor()) // inner and periodic intersections
        stencils_[dd][index] = source_values_[grid_view_.indexSet().index(intersection.outside())];
      else if (!intersection.neighbor() && !intersection.boundary()) // processor boundary
        return false;
      else
        DUNE_THROW(Dune::NotImplemented, "This should not happen!");
    } // intersections
    return true;
  } // void fill_stencils(...)

  const GV& grid_view_;
  const std::vector<LocalVectorType>& source_values_;
  const BoundaryValueType& boundary_values_;
  const AnalyticalFluxType& analytical_flux_;
  std::unique_ptr<SlopeType> slope_;
  const XT::Common::Parameter& param_;
  EigenvectorWrapperType eigenvector_wrapper_;
  DomainType x_local_;
  std::vector<RangeType> slopes_;
  StencilsType stencils_;
};

template <class AnalyticalFluxType,
          class BoundaryValueType,
          class GV,
          class EigenvectorWrapperType = internal::EigenvectorWrapper<AnalyticalFluxType>>
class LocalPointwiseLinearReconstructionOperator : public XT::Grid::ElementFunctor<GV>
{
  using ThisType = LocalPointwiseLinearReconstructionOperator;
  using BaseType = XT::Grid::ElementFunctor<GV>;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using EntityType = typename GV::template Codim<0>::Entity;
  using DomainType = typename BoundaryValueType::DomainType;
  using RangeType = typename BoundaryValueType::RangeReturnType;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using SlopeType = SlopeBase<EntityType, EigenvectorWrapperType>;
  using SlopeFunctorType = LinearSlopeElementFunctor<AnalyticalFluxType, BoundaryValueType, GV, EigenvectorWrapperType>;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, RangeFieldType>;

public:
  explicit LocalPointwiseLinearReconstructionOperator(ReconstructedFunctionType& reconstructed_function,
                                                      const GV& grid_view,
                                                      const std::vector<LocalVectorType>& source_values,
                                                      const BoundaryValueType& boundary_values,
                                                      const AnalyticalFluxType& analytical_flux,
                                                      const SlopeType& slope,
                                                      const XT::Common::Parameter& param,
                                                      const bool flux_is_affine = false)
    : slope_functor_(std::make_unique<SlopeFunctorType>(
          grid_view, source_values, boundary_values, analytical_flux, slope, param, flux_is_affine))
    , reconstructed_function_(reconstructed_function)
  {}

  LocalPointwiseLinearReconstructionOperator(const ThisType& other)
    : BaseType(other)
    , slope_functor_(other.slope_functor_->copy_derived())
    , reconstructed_function_(other.reconstructed_function_)
  {}

  XT::Grid::ElementFunctor<GV>* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    slope_functor_->apply_local(entity);
    // reconstructed function is f(x) = u_entity + slope_matrix * (x - (0.5, 0.5, 0.5, ...))
    auto& local_reconstructed_values = reconstructed_function_.local_values(entity);
    const RangeType u_entity = slope_functor_->u_entity();
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t ii = 0; ii < 2; ++ii) {
        DomainType coord(0.5);
        coord[dd] = ii;
        local_reconstructed_values[coord] = u_entity + slope_functor_->slopes()[dd] * (ii - 0.5);
      } // ii
    } // dd
  } // void apply_local(...)

private:
  std::unique_ptr<SlopeFunctorType> slope_functor_;
  ReconstructedFunctionType& reconstructed_function_;
}; // class LocalPointwiseLinearReconstructionOperator


template <class AnalyticalFluxType,
          class BoundaryValueType,
          class GV,
          class TargetSpaceType,
          class TargetVectorType = typename XT::LA::Container<typename AnalyticalFluxType::R>::VectorType,
          class EigenvectorWrapperType = internal::EigenvectorWrapper<AnalyticalFluxType>>
class LocalLinearReconstructionOperator : public XT::Grid::ElementFunctor<GV>
{
  using BaseType = XT::Grid::ElementFunctor<GV>;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using EntityType = typename GV::template Codim<0>::Entity;
  using RangeType = typename BoundaryValueType::RangeReturnType;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using SlopeType = SlopeBase<EntityType, EigenvectorWrapperType>;
  using TargetType = DiscreteFunction<TargetVectorType, GV, dimRange, 1, RangeFieldType>;
  using TargetBasisType = typename TargetType::SpaceType::GlobalBasisType::LocalizedType;
  using LocalDofVectorType = typename TargetType::DofVectorType::LocalDofVectorType;
  using SlopeFunctorType = LinearSlopeElementFunctor<AnalyticalFluxType, BoundaryValueType, GV, EigenvectorWrapperType>;

public:
  explicit LocalLinearReconstructionOperator(const std::vector<LocalVectorType>& source_values,
                                             TargetSpaceType target_space,
                                             TargetVectorType& target_vector,
                                             const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const SlopeType& slope,
                                             const XT::Common::Parameter& param,
                                             const bool flux_is_affine = false)
    : slope_functor_(std::make_unique<SlopeFunctorType>(
          target_space.grid_view(), source_values, boundary_values, analytical_flux, slope, param, flux_is_affine))
    , target_space_(target_space)
    , target_vector_(target_vector)
    , target_(target_space_, target_vector_, "range")
    , target_basis_(target_.space().basis().localize())
    , local_dof_vector_(target_.dofs().localize())
  {}

  LocalLinearReconstructionOperator(const LocalLinearReconstructionOperator& other)
    : BaseType(other)
    , slope_functor_(other.slope_functor_->copy_derived())
    , target_space_(other.target_space_)
    , target_vector_(other.target_vector_)
    , target_(target_space_, target_vector_, "range")
    , target_basis_(target_.space().basis().localize())
    , local_dof_vector_(target_.dofs().localize())
  {}

  XT::Grid::ElementFunctor<GV>* copy() override final
  {
    return new LocalLinearReconstructionOperator(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    slope_functor_->apply_local(entity);
    // reconstructed function is f(x) = u_entity + slope_matrix * (x - (0.5, 0.5, 0.5, ...))
    local_dof_vector_.bind(entity);
    target_basis_->bind(entity);
    const RangeType u_entity = slope_functor_->u_entity();
    target_basis_->interpolate(
        [&](const auto& xx) {
          auto ret = u_entity;
          for (size_t dd = 0; dd < dimDomain; ++dd)
            ret += slope_functor_->slopes()[dd] * (xx[dd] - 0.5);
          return ret;
        },
        1,
        local_dof_vector_);
  } // void apply_local(...)

private:
  std::unique_ptr<SlopeFunctorType> slope_functor_;
  TargetSpaceType target_space_;
  TargetVectorType& target_vector_;
  TargetType target_;
  std::unique_ptr<TargetBasisType> target_basis_;
  LocalDofVectorType local_dof_vector_;
}; // class LocalLinearReconstructionOperator

template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class GV,
          class MatrixImp = typename XT::LA::Container<typename AnalyticalFluxImp::R>::MatrixType,
          class EigenvectorWrapperImp = internal::EigenvectorWrapper<
              AnalyticalFluxImp,
              FieldMatrix<typename BoundaryValueImp::R, BoundaryValueImp::r, BoundaryValueImp::r>,
              FieldVector<typename BoundaryValueImp::R, BoundaryValueImp::r>>>
class LinearReconstructionOperator : public OperatorInterface<MatrixImp, GV, BoundaryValueImp::r>
{
  using BaseType = OperatorInterface<MatrixImp, GV, BoundaryValueImp::r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  using EigenvectorWrapperType = EigenvectorWrapperImp;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using MatrixType = typename EigenvectorWrapperType::MatrixType;
  using E = XT::Grid::extract_entity_t<GV>;
  using SlopeType = SlopeBase<E, EigenvectorWrapperType>;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using ReconstructionSpaceType = DiscontinuousLagrangeSpace<GV, dimRange, typename AnalyticalFluxType::R>;

  LinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const SourceSpaceType& source_space,
                               const SlopeType& slope = default_minmod_slope(),
                               const bool flux_is_affine = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , source_space_(source_space)
    , range_space_(source_space_.grid_view(), 1)
    , slope_(slope)
    , flux_is_affine_(flux_is_affine)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override
  {
    // evaluate cell averages
    const auto& grid_view = source_space_.grid_view();
    const auto& index_set = grid_view.indexSet();
    std::vector<LocalVectorType> source_values(index_set.size(0));
    auto source_func = make_discrete_function(source_space(), source);
    const auto local_source = source_func.local_function();
    for (const auto& entity : Dune::elements(grid_view)) {
      local_source->bind(entity);
      const auto entity_index = index_set.index(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator = LocalLinearReconstructionOperator<AnalyticalFluxType,
                                                                           BoundaryValueType,
                                                                           GV,
                                                                           ReconstructionSpaceType,
                                                                           VectorType,
                                                                           EigenvectorWrapperType>(
        source_values, range_space_, range, analytical_flux_, boundary_values_, slope_, param, flux_is_affine_);
    auto walker = XT::Grid::Walker<GV>(grid_view);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  static SlopeType& default_minmod_slope()
  {
    static MinmodSlope<E, EigenvectorWrapperType> minmod_slope_;
    return minmod_slope_;
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SourceSpaceType& source_space_;
  ReconstructionSpaceType range_space_;
  const SlopeType& slope_;
  const bool flux_is_affine_;
}; // class LinearReconstructionOperator<...>


// Does not reconstruct a full first-order DG function, but only stores the reconstructed values at the intersection
// centers. This avoids the interpolation in this operator and the evaluation of the reconstructed function in the
// finite volume operator which are both quite expensive for large dimRange.
template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class GV,
          class VectorType,
          class EigenvectorWrapperImp = internal::EigenvectorWrapper<
              AnalyticalFluxImp,
              FieldMatrix<typename BoundaryValueImp::R, BoundaryValueImp::r, BoundaryValueImp::r>,
              FieldVector<typename BoundaryValueImp::R, BoundaryValueImp::r>>>
class PointwiseLinearReconstructionOperator
{
public:
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  using EigenvectorWrapperType = EigenvectorWrapperImp;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using E = XT::Grid::extract_entity_t<GV>;
  using SlopeType = SlopeBase<E, EigenvectorWrapperType>;
  using R = typename BoundaryValueType::R;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using SpaceType = SpaceInterface<GV, dimRange, 1, R>;
  using ReconstructedFunctionType = DiscreteValuedGridFunction<GV, dimRange, 1, R>;
  using ReconstructedValuesType = std::vector<typename ReconstructedFunctionType::LocalFunctionValuesType>;

  PointwiseLinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                                        const BoundaryValueType& boundary_values,
                                        const SpaceType& space,
                                        const SlopeType& slope = default_minmod_slope(),
                                        const bool flux_is_affine = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , space_(space)
    , slope_(slope)
    , flux_is_affine_(flux_is_affine)
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
    auto source_func = make_discrete_function(space_, source);
    const auto local_source = source_func.local_function();
    for (const auto& entity : Dune::elements(grid_view)) {
      local_source->bind(entity);
      const auto entity_index = index_set.index(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator =
        LocalPointwiseLinearReconstructionOperator<AnalyticalFluxType, BoundaryValueType, GV, EigenvectorWrapperType>(
            range, grid_view, source_values, boundary_values_, analytical_flux_, slope_, param, flux_is_affine_);
    auto walker = XT::Grid::Walker<GV>(grid_view);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  static SlopeType& default_minmod_slope()
  {
    static MinmodSlope<E, EigenvectorWrapperType> minmod_slope_;
    return minmod_slope_;
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SpaceType& space_;
  const SlopeType& slope_;
  const bool flux_is_affine_;
}; // class PointwiseLinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_RECONSTRUCTION_LINEAR_HH
