// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_MUSTA_HH
#define DUNE_GDT_OPERATORS_FV_MUSTA_HH

#include <dune/gdt/local/fluxes/musta.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {

template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp, class Traits>
class AdvectionMustaOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, class LocalizableFunctionImp>
class AdvectionMustaOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  typedef typename Dune::GDT::MustaLocalNumericalCouplingFlux<AnalyticalFluxType, LocalizableFunctionType>
      NumericalCouplingFluxType;
  typedef typename Dune::GDT::
      MustaLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType, LocalizableFunctionType>
          NumericalBoundaryFluxType;
  typedef AdvectionMustaOperator<AnalyticalFluxImp,
                                 BoundaryValueImp,
                                 LocalizableFunctionImp,
                                 AdvectionMustaOperatorTraits>
      derived_type;
}; // class AdvectionMustaOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class LocalizableFunctionImp,
          class Traits =
              internal::AdvectionMustaOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, LocalizableFunctionImp>>
class AdvectionMustaOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DomainType;
  using typename BaseType::OnedQuadratureType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  AdvectionMustaOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const size_t num_stages = 2)
    : BaseType(analytical_flux, boundary_values)
    , dx_(dx)
    , num_stages_(num_stages)
  {
  }

  AdvectionMustaOperator(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const OnedQuadratureType& quadrature_1d,
                         const size_t num_stages = 2)
    : BaseType(analytical_flux, boundary_values, quadrature_1d)
    , dx_(dx)
    , num_stages_(num_stages)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, dx_, num_stages_);
  }

private:
  const LocalizableFunctionType& dx_;
  const size_t num_stages_;
}; // class AdvectionMustaOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_MUSTA_HH
