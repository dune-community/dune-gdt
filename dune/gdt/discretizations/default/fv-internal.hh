// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH

#include <dune/xt/common/memory.hh>

#include <dune/gdt/operators/fv.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod,
          class RealizabilityLimiterType = NonLimitingRealizabilityLimiter<typename AnalyticalFluxType::EntityType>>
struct AdvectionOperatorCreator
{
  typedef AdvectionGodunovOperator<AnalyticalFluxType,
                                   BoundaryValueType,
                                   reconstruction_order,
                                   slope_limiter,
                                   RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& /*dx_function*/,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, quadrature_1d, realizability_limiter, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., godunov, ...>

template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::laxfriedrichs,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                         BoundaryValueType,
                                         ConstantFunctionType,
                                         reconstruction_order,
                                         slope_limiter,
                                         RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, dx_function, quadrature_1d, realizability_limiter, false, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., laxfriedrichs, ...>

template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::local_laxfriedrichs,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType,
                                         BoundaryValueType,
                                         ConstantFunctionType,
                                         reconstruction_order,
                                         slope_limiter,
                                         RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, dx_function, quadrature_1d, realizability_limiter, true, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., local_laxfriedrichs, ...>

template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::force,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionForceOperator<AnalyticalFluxType,
                                 BoundaryValueType,
                                 ConstantFunctionType,
                                 reconstruction_order,
                                 slope_limiter,
                                 RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, dx_function, quadrature_1d, realizability_limiter, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., force, ...>


template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::musta,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionMustaOperator<AnalyticalFluxType,
                                 BoundaryValueType,
                                 ConstantFunctionType,
                                 reconstruction_order,
                                 slope_limiter,
                                 RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, dx_function, quadrature_1d, realizability_limiter, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., musta, ...>


template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::laxwendroff,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionLaxWendroffOperator<AnalyticalFluxType,
                                       BoundaryValueType,
                                       ConstantFunctionType,
                                       reconstruction_order,
                                       slope_limiter,
                                       RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, dx_function, quadrature_1d, realizability_limiter, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., laxwendroff, ...>

#if 0
template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType,
          size_t reconstruction_order,
          SlopeLimiters slope_limiter,
          class RealizabilityLimiterType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::kinetic,
                                reconstruction_order,
                                slope_limiter,
                                RealizabilityLimiterType>
{
  typedef AdvectionKineticOperator<AnalyticalFluxType,
                                   BoundaryValueType,
                                   reconstruction_order,
                                   slope_limiter,
                                   RealizabilityLimiterType>
      type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function,
                                      const bool is_linear,
                                      const OnedQuadratureType& quadrature_1d = type::default_1d_quadrature(),
                                      const std::shared_ptr<RealizabilityLimiterType>& realizability_limiter = nullptr)
  {
    return XT::Common::make_unique<type>(
        analytical_flux, boundary_values, quadrature_1d, realizability_limiter, is_linear);
  }
}; // struct AdvectionOperatorCreator<..., kinetic, ...>
#endif


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH
