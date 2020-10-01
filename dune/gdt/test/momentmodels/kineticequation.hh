// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH

#if HAVE_DUNE_XT_DATA

#  include <dune/xt/grid/gridprovider.hh>

#  include <dune/xt/functions/generic/flux-function.hh>
#  include <dune/xt/functions/generic/function.hh>

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class KineticEquationInterface
{
  using ThisType = KineticEquationInterface;

public:
  using MomentBasis = MomentBasisImp;
  using DomainFieldType = typename MomentBasis::DomainFieldType;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  static constexpr size_t dimDomain = MomentBasis::dimDomain;
  static constexpr size_t dimRange = MomentBasis::dimRange;
  static constexpr size_t dimRangeCols = MomentBasis::dimRangeCols;
  static constexpr size_t dimFlux = MomentBasis::dimFlux;
  using FluxType = XT::Functions::FluxFunctionInterface<E, dimRange, dimFlux, dimRange, RangeFieldType>;
  using GenericFluxFunctionType = XT::Functions::GenericFluxFunction<E, dimRange, dimFlux, dimRange, RangeFieldType>;
  using InitialValueType = XT::Functions::FunctionInterface<dimFlux, dimRange, 1, RangeFieldType>;
  using GenericFunctionType = XT::Functions::GenericFunction<dimFlux, dimRange, 1, RangeFieldType>;
  using ScalarFunctionType = XT::Functions::FunctionInterface<dimFlux, 1, 1, RangeFieldType>;
  using BoundaryValueType = InitialValueType;
  using MatrixType = typename Dune::DynamicMatrix<RangeFieldType>;
  using DomainType = typename InitialValueType::DomainType;
  using BasisDomainType = typename MomentBasis::DomainType;
  using StateType = typename FluxType::StateType;
  using RangeReturnType = typename InitialValueType::RangeReturnType;
  using DynamicRangeType = Dune::DynamicVector<RangeFieldType>;

  KineticEquationInterface(const MomentBasis& basis_functions)
    : basis_functions_(basis_functions)
  {}

  virtual ~KineticEquationInterface() {}

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

  virtual std::unique_ptr<FluxType> flux() const = 0;

  virtual std::unique_ptr<InitialValueType> initial_values() const = 0;

  virtual std::unique_ptr<BoundaryValueType> boundary_values() const = 0;

  virtual RangeFieldType CFL() const = 0;

  virtual RangeFieldType t_end() const = 0;

  virtual std::unique_ptr<ScalarFunctionType> sigma_a() const = 0;

  virtual std::unique_ptr<ScalarFunctionType> sigma_s() const = 0;

  virtual std::unique_ptr<ScalarFunctionType> Q() const = 0;

  virtual XT::Common::Configuration grid_config() const = 0;

  static std::string static_id()
  {
    return "kineticequationinterface";
  }

  const MomentBasis& basis_functions() const
  {
    return basis_functions_;
  }

protected:
  const MomentBasis& basis_functions_;
}; // class KineticEquationInterface<E, ...>


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_KINETICEQUATION_HH
