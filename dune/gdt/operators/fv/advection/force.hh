// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_FORCE_HH
#define DUNE_GDT_OPERATORS_FV_FORCE_HH

#include <dune/gdt/local/fluxes/force.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp, class Traits>
class AdvectionForceOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp>
class AdvectionForceOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename Dune::GDT::ForceLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      ForceLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, LocalizableFunctionType>
          NumericalBoundaryFluxType;
  typedef AdvectionForceOperator<AnalyticalFluxImp,
                                 BoundaryValueImp,
                                 LocalizableFunctionImp,
                                 AdvectionForceOperatorTraits>
      derived_type;
}; // class AdvectionForceOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          class Traits =
              internal::AdvectionForceOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, LocalizableFunctionImp>>
class AdvectionForceOperator
  : public Dune::GDT::OperatorInterface<Traits>
  , public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::OnedQuadratureType;
  using typename BaseType::RangeFieldType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionForceOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx)
    : BaseType(analytical_flux, boundary_values)
    , dx_(dx)
  {}

  AdvectionForceOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const OnedQuadratureType& quadrature_1d,
                         const RangeFieldType alpha = BoundaryValueImp::dimDomain)
    : BaseType(analytical_flux, boundary_values, quadrature_1d)
    , dx_(dx)
    , alpha_(alpha)
  {}


  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, alpha_);
  }

private:
  const LocalizableFunctionType& dx_;
  const RangeFieldType alpha_;
}; // class AdvectionForceOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_FORCE_HH
