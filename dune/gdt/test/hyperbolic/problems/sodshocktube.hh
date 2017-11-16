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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SODSHOCKTUBE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SODSHOCKTUBE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/lambda/global-flux-function.hh>
#include <dune/xt/functions/lambda/global-function.hh>

#include <dune/gdt/test/instationary-testcase.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {


// see Lora-Clavijo, F.D. et al, Exact solution of the 1D Riemann problem in Newtonian and relativistic hydrodynamics,
// Rev.Mex.Fis. E59 (2013) 1, 28-50 arXiv:1303.3999, http://www.scielo.org.mx/pdf/rmfe/v59n1/v59n1a5.pdf
// pp. 32-34
// Solution here is not in primitive variables, i.e. u = (rho, rho v, E).
template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class ShocktubeSolution
    : public XT::Functions::GlobalFunctionInterface<EntityImp, DomainFieldImp, 1, RangeFieldImp, 3, 1>
{
  typedef XT::Functions::GlobalFunctionInterface<EntityImp, DomainFieldImp, 1, RangeFieldImp, 3, 1> BaseType;

public:
  static const bool is_linear = false;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".shocktubesolution";
  }

  explicit ShocktubeSolution(const DomainType lower_left,
                             const DomainType upper_right,
                             const XT::Common::ParameterType param_type = {"t", 1},
                             const std::string nm = static_id())
    : lower_left_(lower_left)
    , upper_right_(upper_right)
    , param_type_(param_type)
    , name_(nm)
    , p_ast_(0.303130178050646)
    , v_ast_(-5 * std::sqrt(1.4) * (std::pow(p_ast_, 1.0 / 7.0) - 1))
    , rho_3_(std::pow(p_ast_, 5.0 / 7.0))
    , rho_4_(0.125 * (0.04 + 2.4 * p_ast_) / (0.24 + 0.4 * p_ast_))
    , V_head_(-1.0 * std::sqrt(1.4))
    , V_tail_(v_ast_ - std::sqrt(p_ast_ / rho_3_ * 1.4))
    , V_contact_(v_ast_)
    , V_shock_(std::sqrt(1.12 * (2.4 / 0.28 * p_ast_ + 1.0 / 7.0)))
  {
  }

  ShocktubeSolution(const ShocktubeSolution& other) = default;

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return 1; // TODO: should be at least 7
  }

  virtual void evaluate(const DomainType& x, RangeType& ret, const XT::Common::Parameter& mu) const override final
  {
    double t = this->parse_parameter(mu).get("t")[0];
    if (Dune::XT::Common::FloatCmp::eq(t, 0.0)) {
      if (x < 0.5)
        evaluate_region_1(x[0], ret);
      else
        evaluate_region_6(x[0], ret);
    } else {
      const DomainFieldType x_minus_x0_divided_t = (x[0] - 0.5) / t;
      if (x_minus_x0_divided_t < V_head_)
        evaluate_region_1(x_minus_x0_divided_t, ret);
      else if (x_minus_x0_divided_t < V_tail_)
        evaluate_region_2(x_minus_x0_divided_t, ret);
      else if (x_minus_x0_divided_t < V_contact_)
        evaluate_region_3(x_minus_x0_divided_t, ret);
      else if (x_minus_x0_divided_t < V_shock_)
        evaluate_region_4(x_minus_x0_divided_t, ret);
      else
        evaluate_region_6(x_minus_x0_divided_t, ret);
    }
  }

  virtual void
  jacobian(const DomainType& /*x*/, JacobianRangeType& ret, const XT::Common::Parameter& = {}) const override final
  {
    ret = JacobianRangeType(1);
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual const XT::Common::ParameterType& parameter_type() const override
  {
    return param_type_;
  }

private:
  virtual void evaluate_region_1(const DomainFieldType& /*x*/, RangeType& ret) const
  {
    ret[0] = 1;
    ret[1] = 0;
    ret[2] = 2.5;
  }

  virtual void evaluate_region_2(const DomainFieldType& x, RangeType& ret) const
  {
    const RangeFieldType rho_exact = std::pow(5.0 / 6.0 - 1.0 / (6.0 * std::sqrt(1.4)) * x, 5);
    const RangeFieldType p_exact = std::pow(5.0 / 6.0 - 1.0 / (6.0 * std::sqrt(1.4)) * x, 7);
    const RangeFieldType v_exact = 5.0 / 6.0 * (std::sqrt(1.4) + x);
    ret[0] = rho_exact;
    ret[1] = rho_exact * v_exact;
    ret[2] = 2.5 * p_exact + 0.5 * rho_exact * v_exact * v_exact;
  }

  virtual void evaluate_region_3(const DomainFieldType& /*x*/, RangeType& ret) const
  {
    ret[0] = rho_3_;
    ret[1] = rho_3_ * v_ast_;
    ret[2] = 2.5 * p_ast_ + 0.5 * rho_3_ * v_ast_ * v_ast_;
  }

  virtual void evaluate_region_4(const DomainFieldType& /*x*/, RangeType& ret) const
  {
    ret[0] = rho_4_;
    ret[1] = rho_4_ * v_ast_;
    ret[2] = 2.5 * p_ast_ + 0.5 * rho_4_ * v_ast_ * v_ast_;
  }

  virtual void evaluate_region_6(const DomainFieldType& /*x*/, RangeType& ret) const
  {
    ret[0] = 0.125;
    ret[1] = 0;
    ret[2] = 0.25;
  }

  const DomainType lower_left_;
  const DomainType upper_right_;
  const XT::Common::ParameterType param_type_;
  const std::string name_;
  const double p_ast_;
  const double v_ast_;
  const double rho_3_;
  const double rho_4_;
  const DomainType V_head_;
  const DomainType V_tail_;
  const DomainType V_contact_;
  const DomainType V_shock_;
};


namespace Problems {


template <class E, class D, class U>
class ShockTube : public ProblemBase<E, D, 1, U, typename U::RangeFieldType, 3>
{
  typedef ShockTube<E, D, U> ThisType;
  typedef ProblemBase<E, D, 1, U, typename U::RangeFieldType, 3> BaseType;

public:
  static const bool linear = false;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateRangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;

  typedef typename XT::Functions::GlobalLambdaFluxFunction<U, 0, RangeFieldType, dimRange, dimDomain> ActualFluxType;
  typedef typename XT::Functions::AffineFluxFunction<E, D, dimDomain, U, RangeFieldType, dimRange, 1> ActualRhsType;
  typedef XT::Functions::GlobalLambdaFunction<E, D, dimDomain, RangeFieldType, dimRange, 1> ActualBoundaryValueType;
  typedef XT::Functions::CheckerboardFunction<E, D, dimDomain, RangeFieldType, dimRange, 1> ActualInitialValueType;

  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

  using typename BaseType::FluxType;
  using typename BaseType::RhsType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration grid_config;
    grid_config["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    grid_config["lower_left"] = "[0]";
    grid_config["upper_right"] = "[1]";
    grid_config["num_elements"] = "[8]";
    grid_config["overlap_size"] = "[1]";
    return grid_config;
  }

  static XT::Common::Configuration default_boundary_cfg()
  {
    XT::Common::Configuration boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  ShockTube(const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
            const XT::Common::Configuration& boundary_cfg = default_boundary_cfg())
    : BaseType(create_flux(),
               create_rhs(),
               create_initial_values(grid_cfg),
               create_boundary_values(),
               grid_cfg,
               boundary_cfg,
               0.4,
               0.25,
               false)
  {
  }

  static std::string static_id()
  {
    return "ShockTube";
  }

  static FluxType* create_flux()
  {
    return new ActualFluxType(
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename FluxType::RangeType ret;
          ret[0] = u[1];
          ret[1] = 0.8 * std::pow(u[1], 2) / u[0] + 0.4 * u[2];
          ret[2] = 1.4 * u[1] * u[2] / u[0] - 0.2 * std::pow(u[1], 3) / std::pow(u[0], 2);
          return ret;
        },
        {},
        "sodshocktube flux",
        [](const XT::Common::Parameter&) { return 4; },
        [](const DomainType&, const StateRangeType&, const XT::Common::Parameter&) {
          return typename ActualFluxType::PartialXRangeType(0);
        },
        [](const DomainType&, const StateRangeType& u, const XT::Common::Parameter&) {
          typename ActualFluxType::PartialURangeType ret;
          ret[0] = {0., 1., 0.};
          ret[1] = {-0.8 * std::pow(u[1], 2) / std::pow(u[0], 2), 1.6 * u[1] / u[0], 0.4};
          ret[2] = {-1.4 * u[1] * u[2] / std::pow(u[0], 2) + 0.4 * std::pow(u[1], 3) / std::pow(u[0], 3),
                    1.4 * u[2] / u[0] - 0.6 * std::pow(u[1], 2) / std::pow(u[0], 2),
                    1.4 * u[1] / u[0]};
          return ret;
        });
  }

  static RhsType* create_rhs()
  {
    return new ActualRhsType(FieldVector<MatrixType, 1>(MatrixType(0.)));
  } // ... create_rhs(...)

  static InitialValueType* create_initial_values(const XT::Common::Configuration& grid_cfg)
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg["upper_right"]);
    FieldVector<size_t, dimDomain> num_segments(2);
    std::vector<RangeType> initial_vals(2);
    initial_vals[0] = {1, 0, 2.5};
    initial_vals[1] = {0.125, 0, 0.25};
    return new ActualInitialValueType(lower_left, upper_right, num_segments, initial_vals, "initial_values");
  } // ... create_initial_values()

  virtual BoundaryValueType* create_boundary_values()
  {
    return new ActualBoundaryValueType(
        [=](const DomainType& x, const XT::Common::Parameter&) {
          return RangeType{1 - x[0] * 0.875, 0., 2.5 - x[0] * 2.25};
        },
        1);
  } // ... create_boundary_values()
}; // class Shocktube<...>


} // namespace Problems


template <class G, class R = double>
class ShockTubeTestCase
    : public Dune::GDT::Test::
          InstationaryTestCase<G,
                               Problems::ShockTube<
                                   typename G::template Codim<0>::Entity,
                                   typename G::ctype,
                                   typename internal::DiscreteFunctionProvider<G,
                                                                               GDT::SpaceType::product_fv,
                                                                               0,
                                                                               R,
                                                                               3,
                                                                               1,
                                                                               GDT::Backends::gdt>::type>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  typedef
      typename internal::DiscreteFunctionProvider<G, GDT::SpaceType::product_fv, 0, R, 3, 1, GDT::Backends::gdt>::type
          U;
  typedef typename Problems::ShockTube<E, D, U> ProblemType;

private:
  typedef typename Dune::GDT::Test::InstationaryTestCase<G, ProblemType> BaseType;

public:
  static const size_t dimRange = 3;
  static const size_t dimRangeCols = 1;
  using typename BaseType::SolutionType;

  ShockTubeTestCase(const size_t num_refs = 3, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, ProblemType::default_grid_cfg(), num_refs)
    , exact_solution_(std::make_shared<ShocktubeSolution<E, D, R>>(typename Dune::XT::Common::FieldVector<D, d>(0),
                                                                   typename Dune::XT::Common::FieldVector<D, d>(1)))
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return true;
  }

  virtual const std::shared_ptr<const SolutionType> exact_solution() const override final
  {
    return std::shared_ptr<const SolutionType>(exact_solution_);
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================================================+\n"
        << "|+====================================================================================================+|\n"
        << "||  Testcase: Shock Tube                                                                              ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [0, 1]                                                                                   ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                                                  ||\n"
        << "||  flux = [u[1] 0.8*u[1]*u[1]/u[0]+0.4*u[2] 1.4*u[1]*u[2]/u[0]-0.2*u[1]*u[1]*u[1]/(u[0]*u[0])]       ||\n"
        << "||  rhs = 0                                                                                           ||\n"
        << "||  reference solution: semianalytic solution                                                         ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
  std::shared_ptr<const SolutionType> exact_solution_;
}; // class ShockTubeTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SODSHOCKTUBE_HH
