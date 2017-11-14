// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-function.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>

#include "../base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class ImplementationType>
class KineticEquation : public ProblemBase<typename ImplementationType::EntityType,
                                           typename ImplementationType::DomainFieldType,
                                           ImplementationType::dimDomain,
                                           typename ImplementationType::StateType,
                                           typename ImplementationType::RangeFieldType,
                                           ImplementationType::dimRange>
{
protected:
  typedef ProblemBase<typename ImplementationType::EntityType,
                      typename ImplementationType::DomainFieldType,
                      ImplementationType::dimDomain,
                      typename ImplementationType::StateType,
                      typename ImplementationType::RangeFieldType,
                      ImplementationType::dimRange>
      BaseType;


public:
  using typename BaseType::RangeFieldType;
  static constexpr bool linear = ImplementationType::linear;

  KineticEquation(const ImplementationType& implementation)
    : BaseType(implementation.create_flux(),
               implementation.create_rhs(),
               implementation.create_initial_values(),
               implementation.create_boundary_values(),
               implementation.grid_config(),
               implementation.boundary_config(),
               implementation.CFL(),
               implementation.t_end())
  {
  }

  template <class BasisfunctionType, class GridLayerType>
  KineticEquation(const BasisfunctionType& basis_funcs, const GridLayerType& grid_layer)
    : KineticEquation(ImplementationType(basis_funcs, grid_layer))
  {
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    return ImplementationType::default_grid_cfg();
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    return ImplementationType::default_grid_cfg();
  }

  bool has_non_zero_rhs() const
  {
    return true;
  }

  static std::string static_id()
  {
    return ImplementationType::static_id();
  }
}; // class KineticEquation<...>

template <class BasisfunctionImp, class GridLayerImp, class U_, bool linear_ = true>
class KineticEquationImplementation
{
  typedef KineticEquationImplementation<BasisfunctionImp, GridLayerImp, U_, linear_> ThisType;

public:
  static constexpr bool linear = linear_;
  typedef BasisfunctionImp BasisfunctionType;
  typedef GridLayerImp GridLayerType;
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename BasisfunctionType::DomainFieldType DomainFieldType;
  typedef U_ StateType;
  typedef typename BasisfunctionType::RangeFieldType RangeFieldType;
  static const size_t dimDomain = BasisfunctionType::dimFlux;
  static const size_t dimRange = BasisfunctionType::dimRange;

  typedef XT::Functions::GlobalLambdaFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>
      GlobalLambdaFunctionType;
  typedef XT::Functions::GlobalLambdaFluxFunction<U_, 0, RangeFieldType, dimRange, 1> GlobalLambdaFluxFunctionType;

  typedef typename KineticEquation<ThisType>::FluxType FluxType;
  typedef typename KineticEquation<ThisType>::RhsType RhsType;
  typedef typename KineticEquation<ThisType>::InitialValueType InitialValueType;
  typedef typename KineticEquation<ThisType>::BoundaryValueType BoundaryValueType;

  typedef typename XT::Functions::
      AffineFluxFunction<EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange, 1>
          RhsAffineFunctionType;
  typedef typename XT::Functions::
      AffineFluxFunction<EntityType, DomainFieldType, dimDomain, U_, RangeFieldType, dimRange, dimDomain>
          ActualFluxType;
  typedef XT::Functions::
      CheckerboardFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1, RhsAffineFunctionType>
          ActualRhsType;
  typedef XT::Functions::CheckerboardFunction<EntityType,
                                              DomainFieldType,
                                              dimDomain,
                                              RangeFieldType,
                                              dimRange,
                                              1,
                                              GlobalLambdaFunctionType>
      ActualInitialValueType;
  typedef GlobalLambdaFunctionType ActualBoundaryValueType;

  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename RhsAffineFunctionType::DomainType DomainType;
  typedef typename RhsAffineFunctionType::RangeType RangeType;

  KineticEquationImplementation(const BasisfunctionType& basis_functions, const GridLayerType& grid_layer)
    : basis_functions_(basis_functions)
    , grid_layer_(grid_layer)
  {
  }

  virtual ~KineticEquationImplementation()
  {
  }

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[1.0]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  virtual FluxType* create_flux() const = 0;

  virtual RhsType* create_rhs() const = 0;

  virtual InitialValueType* create_initial_values() const = 0;

  virtual BoundaryValueType* create_boundary_values() const = 0;

  virtual RangeFieldType CFL() const = 0;

  virtual RangeFieldType t_end() const = 0;

  virtual XT::Common::Configuration grid_config() const = 0;

  virtual XT::Common::Configuration boundary_config() const = 0;

  static std::string static_id()
  {
    return "kineticequationimplementation";
  }

protected:
  const BasisfunctionType& basis_functions_;
  const GridLayerType& grid_layer_;
}; // class KineticEquationImplementation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
