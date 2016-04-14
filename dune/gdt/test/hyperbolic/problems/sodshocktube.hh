// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SODSHOCKTUBE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SODSHOCKTUBE_HH

#include <memory>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class ShocktubeSolutionAtSpecificTime
    : public DS::GlobalFunctionInterface<EntityImp, DomainFieldImp, 1, RangeFieldImp, 3, 1>
{
  typedef DS::GlobalFunctionInterface<EntityImp, DomainFieldImp, 1, RangeFieldImp, 3, 1> BaseType;
  typedef ShocktubeSolutionAtSpecificTime<EntityImp, DomainFieldImp, RangeFieldImp> ThisType;

public:
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
    return BaseType::static_id() + ".shocktubesolutionatspecifictime";
  }

  explicit ShocktubeSolutionAtSpecificTime(const double t, const DomainType lower_left, const DomainType upper_right,
                                           const std::string nm = static_id())
    : t_(t)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
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

  ShocktubeSolutionAtSpecificTime(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".shocktubesolutionatspecifictime";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    if (DSC::FloatCmp::eq(t_, 0.0)) {
      if (x < 0.5)
        evaluate_region_1(x[0], ret);
      else
        evaluate_region_6(x[0], ret);
    } else {
      const DomainFieldType x_minus_x0_divided_t = (x[0] - 0.5) / t_;
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

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const override final
  {
    ret = JacobianRangeType(1);
  }

  virtual std::string name() const override final
  {
    return name_;
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
    const RangeFieldType p_exact   = std::pow(5.0 / 6.0 - 1.0 / (6.0 * std::sqrt(1.4)) * x, 7);
    const RangeFieldType v_exact   = 5.0 / 6.0 * (std::sqrt(1.4) + x);
    ret[0]                         = rho_exact;
    ret[1]                         = rho_exact * v_exact;
    ret[2]                         = 2.5 * p_exact + 0.5 * rho_exact * v_exact * v_exact;
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

  const double t_;
  const DomainType lower_left_;
  const DomainType upper_right_;
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

// see Lora-Clavijo, F.D. et al, Exact solution of the 1D Riemann problem in Newtonian and relativistic hydrodynamics,
// Rev.Mex.Fis. E59 (2013) 1, 28-50 arXiv:1303.3999, http://www.scielo.org.mx/pdf/rmfe/v59n1/v59n1a5.pdf
// pp. 32-34
// Solution here is not in primitive variables, i.e. u = (rho, rho v, E).
template <class EntityType, class DomainFieldType, class RangeFieldType>
class ShocktubeSolution
    : public DS::TimeDependentFunctionInterface<
          typename DS::LocalizableFunctionInterface<EntityType, DomainFieldType, 1, RangeFieldType, 3, 1>, double>
{
  typedef typename DS::TimeDependentFunctionInterface<
      typename DS::LocalizableFunctionInterface<EntityType, DomainFieldType, 1, RangeFieldType, 3, 1>, double> BaseType;
  using typename BaseType::TimeIndependentFunctionType;
  typedef ShocktubeSolutionAtSpecificTime<EntityType, DomainFieldType, RangeFieldType> SolutionAtSpecificTimeType;

  typedef typename TimeIndependentFunctionType::DomainType DomainType;

public:
  ShocktubeSolution(const DomainType lower_left, const DomainType upper_right)
    : lower_left_(lower_left)
    , upper_right_(upper_right)
  {
  }

  virtual std::unique_ptr<TimeIndependentFunctionType> evaluate_at_time(const double t) const
  {
    return DSC::make_unique<SolutionAtSpecificTimeType>(t, lower_left_, upper_right_);
  }

  virtual std::string type() const
  {
    return "gdt.shocktubesolution";
  }

  virtual std::string name() const
  {
    return "gdt.shocktubesolution";
  }

private:
  const DomainType lower_left_;
  const DomainType upper_right_;
};


namespace Problems {


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class ShockTube : public Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef ShockTube<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::DefaultRHSType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".shocktube";
  }

  std::string type() const override
  {
    return BaseType::type() + ".shocktube";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0]";
    grid_config["upper_right"]  = "[1]";
    grid_config["num_elements"] = "[8]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["variable"]   = "u";
    flux_config["expression"] = "[u[1] 0.8*u[1]*u[1]/u[0]+0.4*u[2] 1.4*u[1]*u[2]/u[0]-0.2*u[1]*u[1]*u[1]/(u[0]*u[0])]";
    flux_config["order"]      = "4";
    flux_config["gradient.0"] = "[0 1 0; -0.8*u[1]*u[1]/(u[0]*u[0]) 1.6*u[1]/u[0] 0.4; "
                                "-1.4*u[1]*u[2]/(u[0]*u[0])+0.4*u[1]*u[1]*u[1]/(u[0]*u[0]*u[0]) "
                                "1.4*u[2]/u[0]-0.6*u[1]*u[1]/(u[0]*u[0]) 1.4*u[1]/u[0]]";
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"]   = "[0.0]";
    rhs_config["upper_right"]  = "[1.0]";
    rhs_config["num_elements"] = "[1]";
    rhs_config["variable"]     = "u";
    rhs_config["values.0"]     = "[0 0 0]";
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"]   = "[0.0]";
    initial_value_config["upper_right"]  = "[1]";
    initial_value_config["num_elements"] = "[2]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = "[1 0 2.5]";
    initial_value_config["values.1"]     = "[0.125 0 0.25]";
    initial_value_config["order"] = "0";
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config    = DefaultBoundaryValueType::default_config();
    boundary_value_config["type"]       = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = "[1-x[0]*0.875 0 2.5-x[0]*2.25]";
    boundary_value_config["order"] = "1";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  ShockTube(const std::shared_ptr<const FluxType> flux, const std::shared_ptr<const RHSType> rhs,
            const std::shared_ptr<const InitialValueType> initial_values, const ConfigType& grid_config,
            const ConfigType& boundary_info, const std::shared_ptr<const BoundaryValueType> boundary_values)
    : BaseType(flux, rhs, initial_values, grid_config, boundary_info, boundary_values)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 0.25;
  }

  virtual bool is_linear() const override
  {
    return false;
  }
};

} // namespace Problems

template <class G, class R = double>
class ShockTubeTestCase
    : public Dune::GDT::Tests::NonStationaryTestCase<G, Problems::ShockTube<typename G::template Codim<0>::Entity,
                                                                            typename G::ctype, G::dimension, R, 3>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d            = G::dimension;
  static const size_t dimRange     = 3;
  static const size_t dimRangeCols = 1;
  typedef typename Problems::ShockTube<E, D, d, R, dimRange> ProblemType;

private:
  typedef typename Dune::GDT::Tests::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  ShockTubeTestCase(const size_t num_refs = 3, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, Stuff::Grid::Providers::Cube<G>::create(ProblemType::default_grid_config())->grid_ptr(),
               num_refs)
    , problem_(*(ProblemType::create(ProblemType::default_config())))
    , exact_solution_(std::make_shared<ShocktubeSolution<E, D, R>>(typename DSC::FieldVector<D, d>(0),
                                                                   typename DSC::FieldVector<D, d>(1)))
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
        << "||  time = [0, " + DSC::toString(BaseType::t_end())
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
