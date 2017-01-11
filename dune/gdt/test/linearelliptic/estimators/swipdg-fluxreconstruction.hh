// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_ESTIMATORS_SWIPDG_FLUXRECONSTRUCTION_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_ESTIMATORS_SWIPDG_FLUXRECONSTRUCTION_HH

#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <dune/xt/functions/ESV2007.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/local/integrands/ESV2007.hh>
#include <dune/gdt/operators/fluxreconstruction.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

namespace Dune {
namespace GDT {
namespace LinearElliptic {
namespace SwipdgFluxreconstrutionEstimators {


static const size_t over_integrate = 2;


static inline std::string
local_nonconformity_ESV2007_id()
{
  return "eta_NC_ESV2007";
}


static inline std::string
local_residual_ESV2007_id()
{
  return "eta_R_ESV2007";
}


static inline std::string
local_diffusive_flux_ESV2007_id()
{
  return "eta_DF_ESV2007";
}


static inline std::string
ESV2007_id()
{
  return "eta_ESV2007";
}


static inline std::string
ESV2007_alternative_summation_id()
{
  return ESV2007_id() + "_alternative_summation";
}


/**
 *  \brief computes the local nonconformity estimator as defined in ESV2007
 */
template <class SpaceType, class VectorType, class DiffusionFactorType, class DiffusionTensorType,
          class GridViewType = typename SpaceType::GridViewType>
class LocalNonconformityESV2007
    : public XT::Grid::Functor::Codim0Return<GridViewType, typename SpaceType::RangeFieldType>
{
  typedef XT::Grid::Functor::Codim0Return<GridViewType, typename SpaceType::RangeFieldType> BaseType;
  typedef LocalNonconformityESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      ThisType;
  typedef ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;
  typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;
  typedef LocalVolumeIntegralOperator<LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType>>
      LocalOperatorType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::ReturnType;

  static std::string id()
  {
    return local_nonconformity_ESV2007_id();
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(), space, vector, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    ThisType estimator(grid_view, space, vector, diffusion_factor, diffusion_tensor, over_int);
    XT::Grid::Walker<GridViewType> grid_walker(grid_view);
    grid_walker.append(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result());
  } // ... estimate(...)

  LocalNonconformityESV2007(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                            const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                            const size_t over_int = over_integrate)
    : grid_view_(grid_view)
    , discrete_solution_(space, vector)
    , oswald_interpolation_(space)
    , difference_(discrete_solution_ - oswald_interpolation_)
    , local_operator_(over_int, diffusion_factor, diffusion_tensor)
    , prepared_(false)
    , result_(0.0)
  {
  }

  virtual void prepare() override final
  {
    if (prepared_)
      return;
    const OswaldInterpolationOperator<GridViewType> oswald_interpolation_operator(grid_view_);
    oswald_interpolation_operator.apply(discrete_solution_, oswald_interpolation_);
    result_   = 0.0;
    prepared_ = true;
  } // ... prepare(...)

  virtual ReturnType compute_locally(const EntityType& entity) override final
  {
    const auto local_difference = difference_.local_function(entity);
    const auto result           = local_operator_.apply2(*local_difference, *local_difference);
    assert(result.rows() >= 1);
    assert(result.cols() >= 1);
    return result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) override final
  {
    result_ += compute_locally(entity);
  }

  virtual ReturnType result() const override
  {
    return result_;
  }

private:
  const GridViewType& grid_view_;
  const ConstDiscreteFunctionType discrete_solution_;
  DiscreteFunctionType oswald_interpolation_;
  const DifferenceType difference_;
  const LocalOperatorType local_operator_;
  bool prepared_;
  ReturnType result_;
}; // class LocalNonconformityESV2007


/**
 *  \brief computes the local residual estimator as defined in ESV2007
 */
template <class SpaceType, class VectorType, class ForceType, class DiffusionFactorType, class DiffusionTensorType,
          class GridViewType = typename SpaceType::GridViewType>
class LocalResidualESV2007 : public XT::Grid::Functor::Codim0Return<GridViewType, typename ForceType::RangeFieldType>
{
  typedef XT::Grid::Functor::Codim0Return<GridViewType, typename ForceType::RangeFieldType> BaseType;
  typedef LocalResidualESV2007<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      ThisType;
  typedef typename ForceType::RangeFieldType RangeFieldType;
  typedef ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;
  typedef DunePdelabRtSpaceWrapper<GridViewType, 0, RangeFieldType, SpaceType::dimDomain> RTN0SpaceType;
  typedef DiscreteFunction<RTN0SpaceType, VectorType> DiffusiveFluxType;
  typedef typename DiffusiveFluxType::DivergenceType DivergenceType;
  typedef typename DivergenceType::DifferenceType DifferenceType;
  typedef typename XT::Functions::ESV2007::CutoffFunction<DiffusionFactorType, DiffusionTensorType> CutoffFunctionType;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<CutoffFunctionType>> LocalOperatorType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::ReturnType;

  static std::string id()
  {
    return local_residual_ESV2007_id();
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector, const ForceType& force,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(), space, vector, force, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const ForceType& force, const DiffusionFactorType& diffusion_factor,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    ThisType estimator(grid_view, space, vector, force, diffusion_factor, diffusion_tensor, over_int);
    XT::Grid::Walker<GridViewType> grid_walker(grid_view);
    grid_walker.append(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result());
  } // ... estimate(...)

  LocalResidualESV2007(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                       const ForceType& force, const DiffusionFactorType& diffusion_factor,
                       const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
    : grid_view_(grid_view)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , over_integrate_(over_int)
    , discrete_solution_(space, vector)
    , rtn0_space_(grid_view_)
    , diffusive_flux_(rtn0_space_)
    , divergence_(diffusive_flux_.divergence())
    , difference_(force - divergence_)
    , cutoff_function_(diffusion_factor_, diffusion_tensor_)
    , local_operator_(over_integrate_, cutoff_function_)
    , prepared_(false)
    , result_(0.0)
  {
  }

  virtual void prepare() override final
  {
    if (prepared_)
      return;
    const DiffusiveFluxReconstructionOperator<GridViewType, DiffusionFactorType, DiffusionTensorType>
        diffusive_flux_reconstruction(grid_view_, diffusion_factor_, diffusion_tensor_, over_integrate_);
    diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
    result_   = 0.0;
    prepared_ = true;
  } // ... prepare(...)

  virtual ReturnType compute_locally(const EntityType& entity) override final
  {
    const auto local_difference = difference_.local_function(entity);
    auto result                 = local_operator_.apply2(*local_difference, *local_difference);
    assert(result.rows() >= 1);
    assert(result.cols() >= 1);
    return result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) override final
  {
    result_ += compute_locally(entity);
  }

  virtual ReturnType result() const override final
  {
    return result_;
  }

private:
  const GridViewType grid_view_;
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const size_t over_integrate_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  DiffusiveFluxType diffusive_flux_;
  const DivergenceType divergence_;
  const DifferenceType difference_;
  const CutoffFunctionType cutoff_function_;
  const LocalOperatorType local_operator_;
  bool prepared_;
  RangeFieldType result_;
}; // class LocalResidualESV2007


/**
 *  \brief computes the local diffusive flux estimator as defined in ESV2007
 */
template <class SpaceType, class VectorType, class DiffusionFactorType, class DiffusionTensorType,
          class GridViewType = typename SpaceType::GridViewType>
class LocalDiffusiveFluxESV2007
    : public XT::Grid::Functor::Codim0Return<typename SpaceType::GridViewType, typename SpaceType::RangeFieldType>
{
  typedef LocalDiffusiveFluxESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      ThisType;
  typedef XT::Grid::Functor::Codim0Return<typename SpaceType::GridViewType, typename SpaceType::RangeFieldType>
      BaseType;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  typedef ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;
  typedef DunePdelabRtSpaceWrapper<GridViewType, 0, RangeFieldType, SpaceType::dimDomain> RTN0SpaceType;
  typedef DiscreteFunction<RTN0SpaceType, VectorType> RTN0DiscreteFunctionType;
  typedef LocalVolumeIntegralOperator<LocalDiffusiveFluxEstimateESV2007Integrand<DiffusionFactorType,
                                                                                 RTN0DiscreteFunctionType,
                                                                                 DiffusionTensorType>>
      LocalOperatorType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::ReturnType;

  static std::string id()
  {
    return local_diffusive_flux_ESV2007_id();
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(), space, vector, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(grid_view, space, vector, diffusion_factor, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(),
                    space,
                    vector,
                    diffusion_factor_norm,
                    diffusion_factor_reconstruction,
                    diffusion_tensor,
                    over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    ThisType estimator(
        grid_view, space, vector, diffusion_factor_norm, diffusion_factor_reconstruction, diffusion_tensor, over_int);
    XT::Grid::Walker<GridViewType> grid_walker(grid_view);
    grid_walker.append(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result());
  } // ... estimate(...)

  LocalDiffusiveFluxESV2007(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                            const DiffusionFactorType& diffusion_factor_norm,
                            const DiffusionFactorType& diffusion_factor_reconstruction,
                            const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
    : grid_view_(grid_view)
    , diffusion_factor_reconstruction_(diffusion_factor_reconstruction)
    , diffusion_tensor_(diffusion_tensor)
    , discrete_solution_(space, vector)
    , rtn0_space_(grid_view_)
    , diffusive_flux_(rtn0_space_)
    , over_int_(over_int)
    , local_operator_(over_int_, diffusion_factor_norm, diffusion_tensor_, diffusive_flux_)
    , prepared_(false)
    , result_(0.0)
  {
  }

  virtual void prepare() override final
  {
    if (prepared_)
      return;
    const DiffusiveFluxReconstructionOperator<GridViewType, DiffusionFactorType, DiffusionTensorType>
        diffusive_flux_reconstruction(grid_view_, diffusion_factor_reconstruction_, diffusion_tensor_, over_int_);
    diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
    result_   = 0.0;
    prepared_ = true;
  } // ... prepare(...)

  virtual ReturnType compute_locally(const EntityType& entity) override final
  {
    const auto local_discrete_solution = discrete_solution_.local_function(entity);
    auto result                        = local_operator_.apply2(*local_discrete_solution, *local_discrete_solution);
    assert(result.rows() >= 1);
    assert(result.cols() >= 1);
    return result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) override final
  {
    result_ += compute_locally(entity);
  }

  virtual ReturnType result() const override final
  {
    return result_;
  }

private:
  const GridViewType& grid_view_;
  const DiffusionFactorType& diffusion_factor_reconstruction_;
  const DiffusionTensorType& diffusion_tensor_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const size_t over_int_;
  const LocalOperatorType local_operator_;
  bool prepared_;
  RangeFieldType result_;
}; // class LocalDiffusiveFluxESV2007


template <class SpaceType, class VectorType, class ForceType, class DiffusionFactorType, class DiffusionTensorType,
          class GridViewType = typename SpaceType::GridViewType>
class ESV2007 : public XT::Grid::Functor::Codim0Return<GridViewType, typename SpaceType::RangeFieldType>
{
  typedef XT::Grid::Functor::Codim0Return<GridViewType, typename SpaceType::RangeFieldType> BaseType;
  typedef ESV2007<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType, GridViewType> ThisType;

  typedef LocalNonconformityESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      LocalNonconformityEstimator;
  typedef LocalResidualESV2007<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      LocalResidualEstimator;
  typedef LocalDiffusiveFluxESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
      LocalDiffusiveFluxEstimator;

public:
  using typename BaseType::EntityType;
  using typename BaseType::ReturnType;

  static std::string id()
  {
    return ESV2007_id();
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector, const ForceType& force,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(), space, vector, force, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector, const ForceType& force,
                             const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(),
                    space,
                    vector,
                    force,
                    diffusion_factor_norm,
                    diffusion_factor_reconstruction,
                    diffusion_tensor,
                    over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const ForceType& force, const DiffusionFactorType& diffusion_factor,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    return estimate(grid_view, space, vector, force, diffusion_factor, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const ForceType& force, const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    ThisType estimator(grid_view,
                       space,
                       vector,
                       force,
                       diffusion_factor_norm,
                       diffusion_factor_reconstruction,
                       diffusion_tensor,
                       over_int);
    XT::Grid::Walker<GridViewType> grid_walker(grid_view);
    grid_walker.append(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result());
  } // ... estimate(...)

  static XT::LA::CommonDenseVector<ReturnType>
  estimate_local(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                 const ForceType& force, const DiffusionFactorType& diffusion_factor_norm,
                 const DiffusionFactorType& diffusion_factor_reconstruction,
                 const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    XT::LA::CommonDenseVector<ReturnType> local_indicators(boost::numeric_cast<size_t>(grid_view.indexSet().size(0)),
                                                           0.0);
    ReturnType eta_squared = 0.0;
    ThisType estimator(grid_view,
                       space,
                       vector,
                       force,
                       diffusion_factor_norm,
                       diffusion_factor_reconstruction,
                       diffusion_tensor,
                       over_int);
    const auto entity_it_end = grid_view.template end<0>();
    for (auto entity_it = grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity       = *entity_it;
      const auto index         = grid_view.indexSet().index(entity);
      const auto eta_t_squared = estimator.compute_locally(entity);
      local_indicators[index]  = eta_t_squared;
      eta_squared += eta_t_squared;
    }
    for (auto& element : local_indicators)
      element /= eta_squared;
    return local_indicators;
  } // ... estimate_local(...)

  ESV2007(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector, const ForceType& force,
          const DiffusionFactorType& diffusion_factor_norm, const DiffusionFactorType& diffusion_factor_reconstruction,
          const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
    : eta_nc_(grid_view, space, vector, diffusion_factor_norm, diffusion_tensor, over_int)
    , eta_r_(grid_view, space, vector, force, diffusion_factor_norm, diffusion_tensor, over_int)
    , eta_df_(grid_view, space, vector, diffusion_factor_norm, diffusion_factor_reconstruction, diffusion_tensor,
              over_int)
    , result_(0.0)
  {
  }

  virtual void prepare() override final
  {
    eta_nc_.prepare();
    eta_r_.prepare();
    eta_df_.prepare();
  }

  virtual ReturnType compute_locally(const EntityType& entity) override
  {
    return eta_nc_.compute_locally(entity)
           + std::pow(std::sqrt(eta_r_.compute_locally(entity)) + std::sqrt(eta_df_.compute_locally(entity)), 2);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    result_ += compute_locally(entity);
  }

  virtual ReturnType result() const override
  {
    return result_;
  }

protected:
  LocalNonconformityEstimator eta_nc_;
  LocalResidualEstimator eta_r_;
  LocalDiffusiveFluxEstimator eta_df_;
  ReturnType result_;
}; // class ESV2007


template <class SpaceType, class VectorType, class ForceType, class DiffusionFactorType, class DiffusionTensorType,
          class GridViewType = typename SpaceType::GridViewType>
class ESV2007AlternativeSummation
    : public ESV2007<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType, GridViewType>
{
  typedef ESV2007AlternativeSummation<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType,
                                      GridViewType>
      ThisType;
  typedef ESV2007<SpaceType, VectorType, ForceType, DiffusionFactorType, DiffusionTensorType, GridViewType> BaseType;

public:
  using typename BaseType::ReturnType;
  using typename BaseType::EntityType;

  static std::string id()
  {
    return ESV2007_alternative_summation_id();
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector, const ForceType& force,
                             const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                             const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(), space, vector, force, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const SpaceType& space, const VectorType& vector, const ForceType& force,
                             const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    return estimate(space.grid_view(),
                    space,
                    vector,
                    force,
                    diffusion_factor_norm,
                    diffusion_factor_reconstruction,
                    diffusion_tensor,
                    over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const ForceType& force, const DiffusionFactorType& diffusion_factor,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    return estimate(grid_view, space, vector, force, diffusion_factor, diffusion_factor, diffusion_tensor, over_int);
  }

  static ReturnType estimate(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                             const ForceType& force, const DiffusionFactorType& diffusion_factor_norm,
                             const DiffusionFactorType& diffusion_factor_reconstruction,
                             const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    ThisType estimator(grid_view,
                       space,
                       vector,
                       force,
                       diffusion_factor_norm,
                       diffusion_factor_reconstruction,
                       diffusion_tensor,
                       over_int);
    XT::Grid::Walker<GridViewType> grid_walker(grid_view);
    grid_walker.append(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result());
  } // ... estimate(...)

  static XT::LA::CommonDenseVector<ReturnType>
  estimate_local(const GridViewType& grid_view, const SpaceType& space, const VectorType& vector,
                 const ForceType& force, const DiffusionFactorType& diffusion_factor_norm,
                 const DiffusionFactorType& diffusion_factor_reconstruction,
                 const DiffusionTensorType& diffusion_tensor, const size_t over_int = over_integrate)
  {
    XT::LA::CommonDenseVector<ReturnType> local_indicators(boost::numeric_cast<size_t>(grid_view.indexSet().size(0)),
                                                           0.0);
    ThisType estimator(grid_view,
                       space,
                       vector,
                       force,
                       diffusion_factor_norm,
                       diffusion_factor_reconstruction,
                       diffusion_tensor,
                       over_int);
    const auto entity_it_end = grid_view.template end<0>();
    for (auto entity_it = grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity      = *entity_it;
      const auto index        = grid_view.indexSet().index(entity);
      local_indicators[index] = 3.0 * estimator.compute_locally(entity);
    }
    for (auto& element : local_indicators)
      element /= estimator.result();
    return local_indicators;
  } // ... estimate_local(...)

  template <class... Args>
  ESV2007AlternativeSummation(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , eta_nc_squared_(0.0)
    , eta_r_squared_(0.0)
    , eta_df_squared_(0.0)
  {
  }

  virtual ReturnType compute_locally(const EntityType& entity) override final
  {
    const auto eta_nc_t_squared = eta_nc_.compute_locally(entity);
    const auto eta_r_t_squared  = eta_r_.compute_locally(entity);
    const auto eta_df_t_squared = eta_df_.compute_locally(entity);
    eta_nc_squared_ += eta_nc_t_squared;
    eta_r_squared_ += eta_r_t_squared;
    eta_df_squared_ += eta_df_t_squared;
    return eta_nc_t_squared + eta_r_t_squared + eta_df_t_squared;
  }

  virtual ReturnType result() const override
  {
    return std::sqrt(eta_nc_squared_) + std::sqrt(eta_r_squared_) + std::sqrt(eta_df_squared_);
  }

private:
  using BaseType::eta_nc_;
  using BaseType::eta_r_;
  using BaseType::eta_df_;
  ReturnType eta_nc_squared_;
  ReturnType eta_r_squared_;
  ReturnType eta_df_squared_;
}; // class ESV2007AlternativeSummation


} // namespace SwipdgFluxreconstrutionEstimators
} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_ESTIMATORS_SWIPDG_FLUXRECONSTRUCTION_HH
