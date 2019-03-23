// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/walker.hh>

#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

#include "slopes.hh"
#include "internal.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxType,
          class BoundaryValueType,
          class GridViewType,
          class TargetSpaceType,
          class TargetVectorType = typename XT::LA::Container<typename AnalyticalFluxType::R>::VectorType,
          class EigenvectorWrapperType = internal::EigenvectorWrapper<AnalyticalFluxType>>
class LocalLinearReconstructionOperator : public XT::Grid::ElementFunctor<GridViewType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using DomainType = typename AnalyticalFluxType::LocalFunctionType::DomainType;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using RangeType = typename BoundaryValueType::RangeReturnType;
  using RangeFieldType = typename BoundaryValueType::RangeFieldType;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using MatrixType = typename EigenvectorWrapperType::MatrixType;
  using StencilType = DynamicVector<LocalVectorType>;
  using StencilsType = std::vector<DynamicVector<LocalVectorType>>;
  using SlopeType = SlopeBase<LocalVectorType, MatrixType>;
  using TargetType = DiscreteFunction<TargetVectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using TargetBasisType = typename TargetType::SpaceType::GlobalBasisType::LocalizedType;
  using LocalDofVectorType = typename TargetType::DofVectorType::LocalDofVectorType;

public:
  explicit LocalLinearReconstructionOperator(const std::vector<LocalVectorType>& source_values,
                                             TargetSpaceType target_space,
                                             TargetVectorType& target_vector,
                                             const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const SlopeType& slope,
                                             const XT::Common::Parameter& param,
                                             XT::Common::PerThreadValue<EigenvectorWrapperType>& eigenvector_wrapper)
    : source_values_(source_values)
    , target_space_(target_space)
    , target_vector_(target_vector)
    , target_(target_space_, target_vector_, "range")
    , target_basis_(target_.space().basis().localize())
    , local_dof_vector_(target_.dofs().localize())
    , grid_view_(target_.space().grid_view())
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , slope_(slope)
    , param_(param)
    , eigenvector_wrapper_(eigenvector_wrapper)
  {}

  LocalLinearReconstructionOperator(const LocalLinearReconstructionOperator& other)
    : source_values_(other.source_values_)
    , target_space_(other.target_space_)
    , target_vector_(other.target_vector_)
    , target_(target_space_, target_vector_, "range")
    , target_basis_(target_.space().basis().localize())
    , local_dof_vector_(target_.dofs().localize())
    , grid_view_(target_.space().grid_view())
    , analytical_flux_(other.analytical_flux_)
    , boundary_values_(other.boundary_values_)
    , slope_(other.slope_)
    , param_(other.param_)
    , eigenvector_wrapper_(other.eigenvector_wrapper_)
  {}

  virtual XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalLinearReconstructionOperator(*this);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    static const size_t stencil_size = 3;
    thread_local StencilsType stencils(dimDomain, StencilType(stencil_size));
    bool valid = fill_stencils(stencils, entity);
    // In a MPI parallel run, if entity is on boundary of overlap, we do not have to reconstruct
    if (!valid)
      return;

    // get eigenvectors of flux
    const auto entity_index = grid_view_.indexSet().index(entity);
    const auto& u_entity = source_values_[entity_index];
    auto& eigvecs = *eigenvector_wrapper_;
    if (analytical_flux_.x_dependent())
      x_local_ = entity.geometry().local(entity.geometry().center());
    eigvecs.compute_eigenvectors(entity, x_local_, u_entity, param_);

    thread_local StencilType stencil_char(stencil_size);
    // reconstructed function is f(x) = u_entity + slope_matrix * (x - (0.5, 0.5, 0.5, ...))
    thread_local std::vector<LocalVectorType> slopes(dimDomain);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // no need to reconstruct in all directions, as we are only regarding the center of the face, which will
      // always have the same value assigned, independent of the slope in the other directions
      for (size_t ii = 0; ii < stencil_size; ++ii) // transform to characteristic variables
        eigvecs.apply_inverse_eigenvectors(dd, stencils[dd][ii], stencil_char[ii]);
      // get the slope in characteristic coordinates
      const auto slope_char = slope_.get(stencils[dd], stencil_char, eigvecs.eigenvectors(dd));
      eigvecs.apply_eigenvectors(dd, slope_char, slopes[dd]);
    } // dd

    local_dof_vector_.bind(entity);
    target_basis_->bind(entity);
    target_basis_->interpolate(
        [&](const auto& xx) {
          auto ret = stencils[0][1];
          for (size_t dd = 0; dd < dimDomain; ++dd)
            ret += slopes[dd] * (xx[dd] - 0.5);
          return ret;
        },
        1,
        local_dof_vector_);
  } // void apply_local(...)

private:
  bool fill_stencils(StencilsType& stencils, const EntityType& entity) const
  {
    const auto entity_index = grid_view_.indexSet().index(entity);
    for (size_t dd = 0; dd < dimDomain; ++dd)
      stencils[dd][1] = source_values_[entity_index];
    for (const auto& intersection : Dune::intersections(grid_view_, entity)) {
      const size_t dd = intersection.indexInInside() / 2;
      const size_t index = (intersection.indexInInside() % 2) * 2;
      if (intersection.boundary() && !intersection.neighbor()) // boundary intersections
        stencils[dd][index] = boundary_values_.evaluate(entity.geometry().local(intersection.geometry().center()));
      else if (intersection.neighbor()) // inner and periodic intersections
        stencils[dd][index] = source_values_[grid_view_.indexSet().index(intersection.outside())];
      else if (!intersection.neighbor() && !intersection.boundary()) // processor boundary
        return false;
      else
        DUNE_THROW(Dune::NotImplemented, "This should not happen!");
    } // intersections
    return true;
  } // void fill_stencils(...)

private:
  const std::vector<LocalVectorType>& source_values_;
  TargetSpaceType target_space_;
  TargetVectorType& target_vector_;
  TargetType target_;
  std::unique_ptr<TargetBasisType> target_basis_;
  LocalDofVectorType local_dof_vector_;
  const GridViewType& grid_view_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SlopeType& slope_;
  const XT::Common::Parameter& param_;
  DomainType x_local_;
  XT::Common::PerThreadValue<EigenvectorWrapperType>& eigenvector_wrapper_;
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
  using GridViewType = GV;
  using EigenvectorWrapperType = EigenvectorWrapperImp;
  using LocalVectorType = typename EigenvectorWrapperType::VectorType;
  using MatrixType = typename EigenvectorWrapperType::MatrixType;
  using SlopeType = SlopeBase<LocalVectorType, MatrixType>;
  static constexpr size_t dimDomain = BoundaryValueType::d;
  static constexpr size_t dimRange = BoundaryValueType::r;
  using ReconstructionSpaceType = DiscontinuousLagrangeSpace<GridViewType, dimRange, typename AnalyticalFluxType::R>;

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
    , eigenvector_wrapper_(std::reference_wrapper<const AnalyticalFluxType>(analytical_flux_), flux_is_affine)
  {}

  virtual bool linear() const override final
  {
    return false;
  }

  virtual const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  virtual const RangeSpaceType& range_space() const override final
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
                                                                           GridViewType,
                                                                           ReconstructionSpaceType,
                                                                           VectorType,
                                                                           EigenvectorWrapperType>(
        source_values, range_space_, range, analytical_flux_, boundary_values_, slope_, param, eigenvector_wrapper_);
    auto walker = XT::Grid::Walker<GridViewType>(grid_view);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  static SlopeType& default_minmod_slope()
  {
    static MinmodSlope<VectorType, MatrixType> minmod_slope_;
    return minmod_slope_;
  }

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const SourceSpaceType& source_space_;
  ReconstructionSpaceType range_space_;
  const SlopeType& slope_;
  mutable XT::Common::PerThreadValue<EigenvectorWrapperType> eigenvector_wrapper_;
}; // class LinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
