// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_BASE_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_BASE_HH

#include <dune/xt/common/memory.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace internal {


template <class GridType,
          GDT::SpaceType space_type,
          int polOrder = (space_type == GDT::SpaceType::fv) ? 0 : 1,
          class RangeFieldType = double,
          size_t dimRange = 1,
          size_t dimRangeCols = 1,
          GDT::Backends space_backend = GDT::default_space_backend,
          XT::LA::Backends container_backend = XT::LA::default_backend,
          XT::Grid::Layers layer_type = XT::Grid::Layers::leaf>
class DiscreteFunctionProvider
{
  typedef typename GDT::
      SpaceProvider<GridType, layer_type, space_type, space_backend, polOrder, RangeFieldType, dimRange, dimRangeCols>
          SpaceProviderType;

public:
  typedef typename SpaceProviderType::type SpaceImp;
  typedef typename SpaceProviderType::GridLayerType GridLayerType;
  typedef typename XT::LA::Container<RangeFieldType, container_backend>::VectorType VectorType;
  typedef DiscreteFunction<SpaceImp, VectorType> type;

  static type create(GridLayerType grd_layer)
  {
    auto space = SpaceProviderType::create(grd_layer);
    return type(space);
  }
}; // class DiscreteFunctionProvider<...>


} // namespace internal


template <class EntityType, class DomainFieldType, int domainDim, class U_, class RangeFieldType, int rangeDim>
class ProblemBase : public ProblemInterface<EntityType, DomainFieldType, domainDim, U_, RangeFieldType, rangeDim>
{
private:
  typedef ProblemInterface<EntityType, DomainFieldType, domainDim, U_, RangeFieldType, rangeDim> BaseType;

public:
  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;

  ProblemBase(const FluxType& _flux,
              const RhsType& _rhs,
              const InitialValueType& _initial_values,
              const BoundaryValueType& _boundary_values,
              XT::Common::Configuration _grid_cfg,
              XT::Common::Configuration _boundary_cfg,
              const RangeFieldType _CFL,
              const RangeFieldType _t_end,
              const bool _has_non_zero_rhs = true)
    : flux_(_flux)
    , rhs_(_rhs)
    , initial_values_(_initial_values)
    , boundary_values_(_boundary_values)
    , grid_cfg_(_grid_cfg)
    , boundary_cfg_(_boundary_cfg)
    , CFL_(_CFL)
    , t_end_(_t_end)
    , has_non_zero_rhs_(_has_non_zero_rhs)
  {
  }

  /**
   * \note Do not manually delete these pointers, they are managed automatically from here on!
   */
  ProblemBase(const FluxType*&& _flux,
              const RhsType*&& _rhs,
              const InitialValueType*&& _initial_values,
              const BoundaryValueType*&& _boundary_values,
              XT::Common::Configuration _grid_cfg,
              XT::Common::Configuration _boundary_cfg,
              const RangeFieldType _CFL,
              const RangeFieldType _t_end,
              const bool _has_non_zero_rhs = false)
    : flux_(std::move(_flux))
    , rhs_(std::move(_rhs))
    , initial_values_(std::move(_initial_values))
    , boundary_values_(std::move(_boundary_values))
    , grid_cfg_(_grid_cfg)
    , boundary_cfg_(_boundary_cfg)
    , CFL_(_CFL)
    , t_end_(_t_end)
    , has_non_zero_rhs_(_has_non_zero_rhs)
  {
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0 0.0 0.0]";
    grid_config["upper_right"] = "[1.0 1.0 1.0]";
    grid_config["num_elements"] = "[8 8 8]";
    grid_config["overlap_size"] = "[1 1 1]";
    return grid_config;
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "periodic";
    return boundary_config;
  }

  virtual const FluxType& flux() const override
  {
    return flux_.access();
  }

  virtual const RhsType& rhs() const override
  {
    return rhs_.access();
  }

  virtual const InitialValueType& initial_values() const override
  {
    return initial_values_.access();
  }

  virtual const BoundaryValueType& boundary_values() const override
  {
    return boundary_values_.access();
  }

  virtual const XT::Common::Configuration& grid_cfg() const override
  {
    return grid_cfg_;
  }

  virtual const XT::Common::Configuration& boundary_cfg() const override
  {
    return boundary_cfg_;
  }

  virtual RangeFieldType CFL() const override
  {
    return CFL_;
  }

  virtual RangeFieldType t_end() const override
  {
    return t_end_;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return has_non_zero_rhs_;
  }

protected:
  const XT::Common::ConstStorageProvider<FluxType> flux_;
  const XT::Common::ConstStorageProvider<RhsType> rhs_;
  const XT::Common::ConstStorageProvider<InitialValueType> initial_values_;
  const XT::Common::ConstStorageProvider<BoundaryValueType> boundary_values_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_cfg_;
  const RangeFieldType CFL_;
  const RangeFieldType t_end_;
  const bool has_non_zero_rhs_;
}; // class ProblemBase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_PROBLEMS_BASE_HH
