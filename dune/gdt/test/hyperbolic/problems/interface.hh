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

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class U_, class RangeFieldImp, size_t rangeDim>
class ProblemInterface
{
  typedef ProblemInterface<EntityImp, DomainFieldImp, domainDim, U_, RangeFieldImp, rangeDim> ThisType;

public:
  typedef EntityImp EntityType;
  typedef DomainFieldImp DomainFieldType;
  static const size_t dimDomain = domainDim;
  typedef RangeFieldImp RangeFieldType;
  static const size_t dimRange = rangeDim;

  typedef XT::Functions::LocalizableFluxFunctionInterface<EntityType,
                                                          DomainFieldType,
                                                          dimDomain,
                                                          U_,
                                                          0,
                                                          RangeFieldType,
                                                          dimRange,
                                                          dimDomain>
      FluxType;
  typedef XT::Functions::
      LocalizableFluxFunctionInterface<EntityType, DomainFieldType, dimDomain, U_, 0, RangeFieldType, dimRange, 1>
          RhsType;
  typedef XT::Functions::
      LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>
          InitialValueType;
  typedef XT::Functions::
      LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>
          BoundaryValueType;

  typedef typename InitialValueType::DomainType DomainType;
  typedef typename InitialValueType::RangeType RangeType;
  typedef U_ StateType;
  typedef typename StateType::RangeType StateRangeType;

  virtual ~ProblemInterface()
  {
  }

  virtual const FluxType& flux() const = 0;

  virtual const RhsType& rhs() const = 0;

  virtual const InitialValueType& initial_values() const = 0;

  virtual const BoundaryValueType& boundary_values() const = 0;

  virtual const XT::Common::Configuration& grid_cfg() const = 0;

  virtual const XT::Common::Configuration& boundary_cfg() const = 0;

  virtual RangeFieldType CFL() const = 0;

  virtual RangeFieldType t_end() const = 0;
}; // ProblemInterface


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH
