// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tim Keil        (2017)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_ASSEMLBER_HH
#define DUNE_GDT_LOCAL_ASSEMLBER_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalOperatorApplicator : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_local_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;

  LocalOperatorApplicator(const GridLayerType& grid_layer,
                          const LocalOperatorType& local_operator,
                          const SourceType& source,
                          RangeType& range,
                          const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where)
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const
  {
    return where_.apply_on(grid_layer, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.apply(source_, *range_.local_discrete_function(entity));
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where_;
}; // class LocalOperatorApplicator


template <class GridViewType, class LocalOperatorType, class SourceType, class RangeType>
class LocalOperatorJacobianAssembler : public XT::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(is_local_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim0Object<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  LocalOperatorJacobianAssembler(const GridViewType& grid_view,
                                 const LocalOperatorType& local_operator,
                                 const SourceType& x,
                                 const SourceType& source,
                                 RangeType& range,
                                 const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , x_(x)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.assemble_jacobian(x_, source_, *range_.local_discrete_function(entity));
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& x_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where_;
}; // class LocalOperatorJacobianAssembler


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalCouplingOperatorApplicator : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_local_coupling_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalCouplingOperatorApplicator(const GridLayerType& grid_layer,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where,
                                  const XT::Common::Parameter& mu = {})
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
    , mu_(mu)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_layer, intersection);
  }

  virtual void
  apply_local(const IntersectionType& intersection, const EntityType& inside_entity, const EntityType& outside_entity)
  {
    local_operator_.apply(source_,
                          intersection,
                          *range_.local_discrete_function(inside_entity),
                          *range_.local_discrete_function(outside_entity),
                          mu_);
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where_;
  const XT::Common::Parameter& mu_;
}; // class LocalCouplingOperatorApplicator


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalBoundaryOperatorApplicator : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_local_boundary_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalBoundaryOperatorApplicator(const GridLayerType& grid_layer,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where)
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_layer, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& /*outside_entity*/)
  {
    local_operator_.apply(source_, intersection, *range_.local_discrete_function(inside_entity));
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where_;
}; // class LocalBoundaryOperatorApplicator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMLBER_HH
