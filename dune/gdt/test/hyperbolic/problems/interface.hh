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

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/operators/fv/boundary.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class U_, class RangeFieldImp, size_t rangeDim>
class ProblemInterface
{
  using ThisType = ProblemInterface<EntityImp, DomainFieldImp, domainDim, U_, RangeFieldImp, rangeDim>;

public:
  using EntityType = EntityImp;
  using DomainFieldType = DomainFieldImp;
  static const size_t dimDomain = domainDim;
  using RangeFieldType = RangeFieldImp;
  static const size_t dimRange = rangeDim;

  using FluxType = XT::Functions::LocalizableFluxFunctionInterface<EntityType,
                                                                   DomainFieldType,
                                                                   dimDomain,
                                                                   U_,
                                                                   0,
                                                                   RangeFieldType,
                                                                   dimRange,
                                                                   dimDomain>;
  using RhsType = XT::Functions::
      LocalizableFluxFunctionInterface<EntityType, DomainFieldType, dimDomain, U_, 0, RangeFieldType, dimRange, 1>;
  using InitialValueType =
      XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using DirichletBoundaryValueType =
      XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;
  using GridLayerType = typename U_::SpaceType::GridLayerType;
  using IntersectionType = typename GridLayerType::Intersection;
  using StateType = U_;
  using StateRangeType = typename StateType::RangeType;
  using RangeType = typename InitialValueType::RangeType;
  using DomainType = typename InitialValueType::DomainType;
  using BoundaryValueType = LocalizableBoundaryValueInterface<EntityType, IntersectionType, RangeType>;
  using SolutionType = DirichletBoundaryValueType;

  virtual ~ProblemInterface() {}

  virtual const FluxType& flux() const = 0;

  virtual const RhsType& rhs() const = 0;

  virtual const InitialValueType& initial_values() const = 0;

  virtual const BoundaryValueType& boundary_values() const = 0;

  virtual const XT::Common::Configuration& grid_cfg() const = 0;

  virtual const XT::Common::Configuration& boundary_cfg() const = 0;

  virtual RangeFieldType CFL() const = 0;

  virtual RangeFieldType t_end() const = 0;

  virtual bool has_non_zero_rhs() const
  {
    return false;
  }
}; // ProblemInterface


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_INTERFACE_HH
