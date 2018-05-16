// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_BASE_HH
#define DUNE_GDT_LOCAL_FLUXES_BASE_HH

#include <tuple>
#include <type_traits>

namespace Dune {
namespace GDT {
namespace internal {


template <class AnalyticalFluxImp>
struct NumericalCouplingFluxTraitsBase
{
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef typename AnalyticalFluxType::EntityType EntityType;
  typedef typename AnalyticalFluxType::DomainFieldType DomainFieldType;
  typedef typename AnalyticalFluxType::DomainType DomainType;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename AnalyticalFluxType::StateRangeType RangeType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
}; // class NumericalCouplingFluxTraitsBase

template <class AnalyticalFluxImp, class BoundaryValueImp>
struct NumericalBoundaryFluxTraitsBase : public NumericalCouplingFluxTraitsBase<AnalyticalFluxImp>
{
  using BoundaryValueType = BoundaryValueImp;
  using LocalBoundaryValueType = typename BoundaryValueType::LocalBoundaryValueInterfaceType;
}; // class NumericalBoundaryFluxTraitsBase


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_BASE_HH
