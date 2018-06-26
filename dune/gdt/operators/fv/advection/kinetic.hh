// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_KINETIC_HH
#define DUNE_GDT_OPERATORS_FV_KINETIC_HH

#include <dune/gdt/local/fluxes/kinetic.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp, class GridLayerImp, class Traits>
class AdvectionKineticOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, class BasisfunctionImp, class GridLayerImp>
class AdvectionKineticOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  typedef typename Dune::GDT::KineticLocalNumericalCouplingFlux<AnalyticalFluxImp, BasisfunctionImp, GridLayerImp>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      KineticLocalNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp, BasisfunctionImp, GridLayerImp>
          NumericalBoundaryFluxType;
  typedef AdvectionKineticOperator<AnalyticalFluxImp,
                                   BoundaryValueImp,
                                   BasisfunctionImp,
                                   GridLayerImp,
                                   AdvectionKineticOperatorTraits>
      derived_type;
}; // class AdvectionKineticOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class BasisfunctionImp,
          class GridLayerImp,
          class Traits = internal::
              AdvectionKineticOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, BasisfunctionImp, GridLayerImp>>
class AdvectionKineticOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::IntersectionQuadratureType;

  AdvectionKineticOperator(const AnalyticalFluxType& analytical_flux,
                           const BoundaryValueType& boundary_values,
                           const BasisfunctionImp& basis_functions,
                           const IntersectionQuadratureType& intersection_quadrature =
                               midpoint_quadrature<DomainFieldType, BaseType::dimDomain>())
    : BaseType(analytical_flux, boundary_values, intersection_quadrature)
    , basis_functions_(basis_functions)
  {
  }

  template <class SourceType, class RangeType>
  void apply(SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, basis_functions_);
  }

private:
  const BasisfunctionImp& basis_functions_;
}; // class AdvectionKineticOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_KINETIC_HH
