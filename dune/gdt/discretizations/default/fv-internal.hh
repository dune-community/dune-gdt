// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH

#include <dune/xt/common/memory.hh>

#include <dune/gdt/operators/fv.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType, NumericalFluxes numerical_flux>
struct AdvectionOperatorCreator
{
  typedef AdvectionGodunovOperator<AnalyticalFluxType, BoundaryValueType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& /*dx_function*/)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values);
  }
}; // struct AdvectionOperatorCreator<..., godunov, ...>

template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::laxfriedrichs>
{
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values, dx_function);
  }
}; // struct AdvectionOperatorCreator<..., laxfriedrichs, ...>

template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::local_laxfriedrichs>
{
  typedef AdvectionLaxFriedrichsOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values, dx_function, true);
  }
}; // struct AdvectionOperatorCreator<..., local_laxfriedrichs, ...>

template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType, NumericalFluxes::force>
{
  typedef AdvectionForceOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values, dx_function);
  }
}; // struct AdvectionOperatorCreator<..., force, ...>


template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType, NumericalFluxes::musta>
{
  typedef AdvectionMustaOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values, dx_function);
  }
}; // struct AdvectionOperatorCreator<..., musta, ...>


template <class AnalyticalFluxType, class BoundaryValueType, class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::laxwendroff>
{
  typedef AdvectionLaxWendroffOperator<AnalyticalFluxType, BoundaryValueType, ConstantFunctionType> type;
  typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values,
                                      const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>(analytical_flux, boundary_values, dx_function);
  }
}; // struct AdvectionOperatorCreator<..., laxwendroff, ...>

#if 0
template <class AnalyticalFluxType,
          class BoundaryValueType,
          class ConstantFunctionType>
struct AdvectionOperatorCreator<AnalyticalFluxType,
                                BoundaryValueType,
                                ConstantFunctionType,
                                NumericalFluxes::kinetic>
{
  typedef AdvectionKineticOperator<AnalyticalFluxType,
                                   BoundaryValueType>
      type; typedef typename type::OnedQuadratureType OnedQuadratureType;

  static std::unique_ptr<type> create(const AnalyticalFluxType& analytical_flux,
                                      const BoundaryValueType& boundary_values, const ConstantFunctionType& dx_function)
  {
    return XT::Common::make_unique<type>( analytical_flux, boundary_values);
  }
}; // struct AdvectionOperatorCreator<..., kinetic, ...>
#endif


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_FV_INTERNAL_HH
