// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_GODUNOV_HH
#define DUNE_GDT_OPERATORS_FV_GODUNOV_HH

#include <dune/gdt/local/fluxes/godunov.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxImp, class BoundaryValueImp, class Traits>
class AdvectionGodunovOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp>
class AdvectionGodunovOperatorTraits : public AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  typedef AdvectionTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  typedef typename Dune::GDT::GodunovLocalNumericalCouplingFlux<AnalyticalFluxImp> NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp, BoundaryValueImp>
      NumericalBoundaryFluxType;
  typedef AdvectionGodunovOperator<AnalyticalFluxImp, BoundaryValueImp, AdvectionGodunovOperatorTraits> derived_type;
}; // class AdvectionGodunovOperatorTraits


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class Traits = internal::AdvectionGodunovOperatorTraits<AnalyticalFluxImp, BoundaryValueImp>>
class AdvectionGodunovOperator : public Dune::GDT::OperatorInterface<Traits>, public AdvectionOperatorBase<Traits>
{
  typedef AdvectionOperatorBase<Traits> BaseType;

public:
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::BoundaryValueType;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  using MatrixType = FieldMatrix<typename AnalyticalFluxType::DomainFieldType, dimRange, dimRange>;
  using VectorType = FieldVector<typename AnalyticalFluxType::RangeFieldType, dimRange>;
  using JacobianWrapperType = internal::GodunovJacobianWrapper<AnalyticalFluxType, MatrixType, VectorType>;

  AdvectionGodunovOperator(const AnalyticalFluxType& analytical_flux, const BoundaryValueType& boundary_values)
    : BaseType(analytical_flux, boundary_values)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    BaseType::apply(source, range, param, jacobian_wrapper_);
  }

private:
  mutable XT::Common::PerThreadValue<JacobianWrapperType> jacobian_wrapper_;
}; // class AdvectionGodunovOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_GODUNOV_HH
