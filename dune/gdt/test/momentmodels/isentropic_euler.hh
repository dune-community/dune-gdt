// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_ISENTROPIC_EULER_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_ISENTROPIC_EULER_HH

#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/functions/generic/flux-function.hh>
#include <dune/xt/functions/generic/function.hh>

namespace Dune {
namespace GDT {


template <class E>
class IsentropicEulerEquations
{
  using ThisType = IsentropicEulerEquations;

public:
  using DomainFieldType = typename E::Geometry::ctype;
  using RangeFieldType = DomainFieldType;
  static const size_t dimDomain = 1;
  static const size_t dimRange = 2;
  static const size_t dimRangeCols = 1;
  using FluxType = XT::Functions::FluxFunctionInterface<E, dimRange, dimDomain, dimRange, RangeFieldType>;
  using GenericFluxFunctionType = XT::Functions::GenericFluxFunction<E, dimRange, dimDomain, dimRange, RangeFieldType>;
  using InitialValueType = XT::Functions::FunctionInterface<dimDomain, dimRange, 1, RangeFieldType>;
  using GenericFunctionType = XT::Functions::GenericFunction<dimDomain, dimRange, 1, RangeFieldType>;
  using ConstantFunctionType = XT::Functions::ConstantFunction<dimDomain, dimRange, 1, RangeFieldType>;
  using ScalarFunctionType = XT::Functions::FunctionInterface<dimDomain, 1, 1, RangeFieldType>;
  using ConstantScalarFunctionType = XT::Functions::ConstantFunction<dimDomain, 1, 1, RangeFieldType>;
  using BoundaryValueType = InitialValueType;
  using MatrixType = typename Dune::DynamicMatrix<RangeFieldType>;
  using DomainType = typename InitialValueType::DomainType;
  using StateType = typename FluxType::StateType;
  using RangeReturnType = typename InitialValueType::RangeReturnType;
  using DynamicFluxRangeType = typename FluxType::LocalFunctionType::DynamicRangeType;
  using DynamicFluxJacobianRangeType = typename FluxType::LocalFunctionType::DynamicJacobianRangeType;

  IsentropicEulerEquations(const XT::Common::Configuration grid_config = default_grid_cfg())
    : grid_config_(grid_config)
  {}

  virtual ~IsentropicEulerEquations() {}

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

  virtual std::unique_ptr<FluxType> flux() const
  {
    const double gamma = 1.;
    auto eval_func =
        [gamma](const DomainType& /*x*/, const StateType& u, DynamicFluxRangeType& ret, const XT::Common::Parameter&) {
          const auto& rho = u[0];
          const auto& v = u[1];
          const auto kappa = std::pow(gamma - 1, 2) / (4. * gamma);
          ret[0][0] = v;
          ret[0][1] = v * v / rho + kappa * std::pow(rho, gamma);
        };
    auto jacobian_func = [gamma](const DomainType& /*x*/,
                                 const StateType& u,
                                 DynamicFluxJacobianRangeType& ret,
                                 const XT::Common::Parameter&) {
      const auto& rho = u[0];
      const auto& v = u[1];
      const auto kappa = std::pow(gamma - 1, 2) / (4. * gamma);
      ret[0][0][0] = 0;
      ret[0][0][1] = 1;
      ret[0][1][0] = -v * v / (rho * rho) + kappa * gamma * std::pow(rho, gamma - 1);
      ret[0][1][1] = 2 * v / rho;
    };
    return std::make_unique<GenericFluxFunctionType>(20,
                                                     GenericFluxFunctionType::default_post_bind_function(),
                                                     eval_func,
                                                     XT::Common::ParameterType{},
                                                     "flux",
                                                     jacobian_func);
  }

  virtual std::unique_ptr<InitialValueType> initial_values() const
  {
    return std::make_unique<GenericFunctionType>(0, [](const DomainType& x, const XT::Common::Parameter&) {
      return x[0] < 0.5 ? RangeReturnType{1., 1.} : RangeReturnType{2., -1.};
    });
  }

  virtual std::unique_ptr<BoundaryValueType> boundary_values() const
  {
    return initial_values();
  }

  virtual RangeFieldType CFL() const
  {
    return 0.49;
  }

  virtual RangeFieldType t_end() const
  {
    return 10.0;
  }

  virtual XT::Common::Configuration grid_config() const
  {
    return grid_config_;
  }

  static std::string static_id()
  {
    return "isentropic_euler_equations";
  }

private:
  const XT::Common::Configuration grid_config_;
}; // class IsentropicEulerEquations<E>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_ISENTROPIC_EULER_HH
