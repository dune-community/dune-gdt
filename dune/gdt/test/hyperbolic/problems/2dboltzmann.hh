// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/xt/common/string.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/**
 * Testcase for the Boltzmann equation in two dimensions,
 * see Section 4.1 of Brunner, Holloway, "Two-dimensional time dependent Riemann solvers for neutron transport", Journal
 * of Computational Physics, Volume 210, Issue 1, 2005
 * http://dx.doi.org/10.1016/j.jcp.2005.04.011
 * */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
class Boltzmann2DLineSource
    : public Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, (momentOrder + 1) * (momentOrder + 2) / 2>
{
  typedef Boltzmann2DLineSource<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> ThisType;
  typedef Default<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, (momentOrder + 1) * (momentOrder + 2) / 2>
      BaseType;

public:
  static const bool linear = true;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DummyEntityType;
  typedef typename Dune::Stuff::Functions::Affine<DummyEntityType, RangeFieldImp, dimRange, RangeFieldImp, dimRange,
                                                  dimDomain>
      FluxAffineFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxAffineFunctionType, EntityImp, DomainFieldImp,
                                                                dimDomain, RangeFieldImp, dimRange, 1>
      DefaultFluxType;
  typedef typename DefaultFluxType::FluxRangeType FluxRangeType;
  typedef typename FluxAffineFunctionType::FieldMatrixType MatrixType;
  using typename BaseType::DefaultInitialValueType;
  typedef typename DS::Functions::Affine<DummyEntityType, RangeFieldImp, dimRange, RangeFieldImp, dimRange, 1>
      RHSAffineFunctionType;
  typedef typename DS::Functions::FunctionCheckerboard<RHSAffineFunctionType, EntityImp, DomainFieldImp, dimDomain,
                                                       RangeFieldImp, dimRange, 1>
      RHSCheckerboardFunctionType;
  typedef typename Dune::GDT::CheckerboardBasedRhsEvaluationFlux<RHSCheckerboardFunctionType, EntityImp, DomainFieldImp,
                                                                 dimDomain, RangeFieldImp, dimRange, 1>
      DefaultRHSType;
  typedef typename DefaultRHSType::RangeType RangeType;
  typedef typename DefaultRHSType::DomainType DomainType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".boltzmann2d";
  }

  std::string type() const
  {
    return BaseType::type() + ".boltzmann2d";
  }

  static std::string short_id()
  {
    return "Boltzmann2DLineSource";
  }

protected:
  class GetData
  {
  public:
    static const int precision = 15; // precision for to_string

    // (\Sigma_s \delta_{l0}\delta{m0} - \Sigma_t) * \psi_l^m
    static void create_rhs_values(ConfigType& rhs_config)
    {
      rhs_config["lower_left"]    = "[0.0 0.0]";
      rhs_config["upper_right"]   = "[1.0 1.0]";
      rhs_config["num_elements"]  = "[1 1]";
      const RangeFieldImp Sigma_s = 0;
      const RangeFieldImp Sigma_t = 0;
      MatrixType S;
      S[pos(0, 0)][pos(0, 0)] = Sigma_s - Sigma_t;
      for (size_t l = 1; l <= momentOrder; ++l)
        for (size_t m = 0; m <= l; ++m)
          S[pos(l, m)][pos(l, m)] = -1.0 * Sigma_t;
      rhs_config["A.0"] = Dune::XT::Common::to_string(S, precision);
      rhs_config["b"]   = Dune::XT::Common::to_string(RangeType(0));
    } // ... create_rhs_values(...)

    static void create_flux_matrices(ConfigType& flux_config)
    {
      // X is the matrix corresponding to the x derivative, Z is corresponding to z derivative
      MatrixType X, Z;
      // l, m as in \psi_l^m
      for (size_t l = 0; l <= momentOrder; ++l) {
        // m = 0 case
        size_t row = pos(l, 0);
        if (l > 1)
          X[row][pos(l - 1, 1)] = E(l - 1, 1);
        if (l > 0)
          Z[row][pos(l - 1, 0)] = A(l - 1, 0);
        if (l < momentOrder) {
          X[row][pos(l + 1, 1)] = -1.0 * F(l + 1, 1);
          Z[row][pos(l + 1, 0)] = B(l + 1, 0);
        }
        // m > 0
        for (size_t m = 1; m <= l; ++m) {
          row = pos(l, m);
          if (l > 0) {
            X[row][pos(l - 1, m - 1)] = -0.5 * C(l - 1, m - 1);
            if (m < l - 1)
              X[row][pos(l - 1, m + 1)] = 0.5 * E(l - 1, m + 1);
            if (m < l)
              Z[row][pos(l - 1, m)] = A(l - 1, m);
          }
          if (l < momentOrder) {
            X[row][pos(l + 1, m - 1)] = 0.5 * D(l + 1, m - 1);
            X[row][pos(l + 1, m + 1)] = -0.5 * F(l + 1, m + 1);
            Z[row][pos(l + 1, m)]     = B(l + 1, m);
          }
        }
      }
      flux_config["A.0"] = Dune::XT::Common::to_string(X, precision);
      flux_config["A.1"] = Dune::XT::Common::to_string(Z, precision);
      flux_config["b"]   = Dune::XT::Common::to_string(FluxRangeType(0));
    } // ... create_flux_matrix()

    // initial value is max(exp(-10*((x-0.5)^2 + (y-0.5)^2)/sigma^2), 10^(-4)) with sigma = 0.02 for \psi_0^0 and 0 else
    // expression does not offer max, so use max(a,b) = 0.5*(a+b+abs(a-b))
    static std::string create_initial_values()
    {
      std::string str = "[";
      str += "0.5*(exp(-10*((x[0]-0.5)^2+(x[1]-0.5)^2)/0.0004)+0.0001+abs(exp(-10*((x[0]-0.5)^2+(x[1]-0.5)^2)/"
             "0.0004)-0.0001))";
      for (size_t ii = 1; ii < dimRange; ++ii)
        str += " 0";
      str += "]";
      return str;
    } // ... create_initial_values()

    // boundary values are 0
    static std::string create_boundary_values()
    {
      std::string str = "[";
      for (size_t ii = 0; ii < dimRange; ++ii) {
        if (ii > 0)
          str += " ";
        str += "0";
      }
      str += "]";
      return str;
    } // ... create_boundary_values()

  protected:
    // calculates position of \psi_l^m in vector.
    // The \psi_l^m are ordered by m first and then by l,
    // i.e. (\psi_0^0, \psi_1^0, \psi_2^0, ..., \psi_N^0, \psi_1^1, ..., \psi_1^1, \psi_2^2, ...)
    // Thus \psi_l^m has position l - m + \sum_{k=0}^(l-1) ((N+1) - k ) = (2N - m + 1)*m/2 + l in the vector.
    static size_t pos(const size_t l = 0, const size_t m = 0)
    {
      return ((2.0 * momentOrder - m + 1.0) * m) / 2.0 + l;
    }

    static RangeFieldImp A(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l - m + 1.0) * (l + m + 1.0)) / ((2.0 * l + 3.0) * (2.0 * l + 1.0)));
    }

    static RangeFieldImp B(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l - m) * (l + m)) / ((2.0 * l + 1.0) * (2.0 * l - 1.0)));
    }

    static RangeFieldImp C(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l + m + 1.0) * (l + m + 2.0)) / ((2.0 * l + 3.0) * (2.0 * l + 1.0)));
    }

    static RangeFieldImp D(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l - m) * (l - m - 1.0)) / ((2.0 * l + 1.0) * (2.0 * l - 1.0)));
    }

    static RangeFieldImp E(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l - m + 1.0) * (l - m + 2.0)) / ((2.0 * l + 3.0) * (2.0 * l + 1.0)));
    }

    static RangeFieldImp F(const size_t l, const size_t m)
    {
      assert(m <= l && l >= 0 && l <= momentOrder);
      return std::sqrt(((l + m) * (l + m - 1.0)) / ((2.0 * l + 1.0) * (2.0 * l - 1.0)));
    }
  }; // class GetData

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0 0.0]";
    grid_config["upper_right"]  = "[1.0 1.0]";
    grid_config["num_elements"] = "[60 60]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static std::unique_ptr<ThisType> create(const ConfigType cfg       = default_config(),
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

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    GetData::create_flux_matrices(flux_config);
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    ;
    GetData::create_rhs_values(rhs_config);
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"]   = "[0.0 0.0]";
    initial_value_config["upper_right"]  = "[1.0 1.0]";
    initial_value_config["num_elements"] = "[1 1]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = GetData::create_initial_values();
    initial_value_config["name"]         = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"]       = BoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = GetData::create_boundary_values();
    boundary_value_config["order"]      = "0";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Boltzmann2DLineSource(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> rhs_in,
                        const std::shared_ptr<const InitialValueType> initial_values_in,
                        const ConfigType& grid_config_in, const ConfigType& boundary_info_in,
                        const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }
}; // ... Boltzmann2DLineSource ...


/**
 * Testcase for the Boltzmann equation in two dimensions,
 * see Section 4.2 of Brunner, Holloway, "Two-dimensional time dependent Riemann solvers for neutron transport", Journal
 * of Computational Physics, Volume 210, Issue 1, 2005
 * http://dx.doi.org/10.1016/j.jcp.2005.04.011
 * */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
class Boltzmann2DCheckerboard
    : public Boltzmann2DLineSource<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder>
{
  typedef Boltzmann2DCheckerboard<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> ThisType;
  typedef Boltzmann2DLineSource<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::RangeType;
  using typename BaseType::DomainType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".boltzmann2dcheckerboard";
  }

  std::string type() const override
  {
    return BaseType::type() + ".boltzmann2dcheckerboard";
  }

  static std::string short_id()
  {
    return "Boltzmann2DCheckerBoard";
  }

protected:
  class GetData : public BaseType::GetData
  {
  public:
    using BaseType::GetData::precision; // precision for to_string

    // rhs is q - (\Sigma_s \delta_{l0}\delta{m0} - \Sigma_t) * \psi_l^m
    // see the Checkerboard test case in http://www.sciencedirect.com/science/article/pii/S0021999105002275?np=y
    // q = 0 except in the center where q = 1. Sigma_s = Sigma_t = 1 in scattering regions, Sigma_s = 0, Sigma_t
    // = 10 in absorbing regions. Center is also a scattering region.
    static bool is_absorbing(const size_t row, const size_t col)
    {
      return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
             || ((row == 3 || row == 5) && (col == 1 || col == 5));
    }

    static void create_rhs_values(ConfigType& rhs_config)
    {
      rhs_config["lower_left"]   = "[0.0 0.0]";
      rhs_config["upper_right"]  = "[7.0 7.0]";
      rhs_config["num_elements"] = "[7 7]";
      MatrixType S;
      RangeFieldImp Sigma_s;
      RangeFieldImp Sigma_t;
      RangeType q;
      for (size_t row = 0; row < 7; ++row) {
        for (size_t col = 0; col < 7; ++col) {
          if (row == 3 && col == 3) { // center
            q *= 0;
            q[0]    = 1;
            Sigma_s = 1;
            Sigma_t = 1;
          } else if (is_absorbing(row, col)) { // absorbing regions
            q *= 0;
            Sigma_s = 0;
            Sigma_t = 10;
          } else { // scattering regions (without center)
            q *= 0;
            Sigma_s = 1;
            Sigma_t = 1;
          }
          S *= 0.0;
          S[pos(0, 0)][pos(0, 0)] = Sigma_s - Sigma_t;
          for (size_t l = 1; l <= momentOrder; ++l)
            for (size_t m = 0; m <= l; ++m)
              S[pos(l, m)][pos(l, m)] = -1.0 * Sigma_t;
          size_t number                                          = 7 * row + col;
          rhs_config["A." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(S, precision);
          rhs_config["b." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(q);
        }
      }
    } // ... create_rhs_values(...)

    // initial value is 0
    static std::string create_initial_values()
    {
      std::string str = "[";
      str += "0";
      for (size_t ii = 1; ii < dimRange; ++ii)
        str += " 0";
      str += "]";
      return str;
    } // ... create_initial_values()

  protected:
    using BaseType::GetData::pos;
  }; // class GetData

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"]         = "provider.cube";
    grid_config["lower_left"]   = "[0.0 0.0]";
    grid_config["upper_right"]  = "[7.0 7.0]";
    grid_config["num_elements"] = "[14 14]";
    return grid_config;
  }

  static ConfigType default_boundary_info_config()
  {
    ConfigType boundary_config;
    boundary_config["type"] = "dirichlet";
    return boundary_config;
  }

  static std::unique_ptr<ThisType> create(const ConfigType cfg       = default_config(),
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

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config = BaseType::default_config().sub("flux");
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    GetData::create_rhs_values(rhs_config);
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"]   = "[0.0 0.0]";
    initial_value_config["upper_right"]  = "[7.0 7.0]";
    initial_value_config["num_elements"] = "[1 1]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = GetData::create_initial_values();
    initial_value_config["name"]         = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config = BaseType::default_config().sub("boundary_values");
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Boltzmann2DCheckerboard(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> rhs_in,
                          const std::shared_ptr<const InitialValueType> initial_values_in,
                          const ConfigType& grid_config_in, const ConfigType& boundary_info_in,
                          const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, rhs_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }

  virtual double CFL() const override
  {
    return 0.4;
  }

  virtual double t_end() const override
  {
    return 3.2;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return true;
  }
}; // ... Boltzmann2DCheckerboard ...


} // namespace Problems


template <class G, class R = double, size_t momentOrder = 1>
class Boltzmann2DCheckerboardTestCase
    : public Dune::GDT::Tests::
          NonStationaryTestCase<G, Problems::Boltzmann2DCheckerboard<typename G::template Codim<0>::Entity,
                                                                     typename G::ctype, G::dimension, R, momentOrder>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 2, "Only implemented for dimension 2.");
  typedef typename Problems::Boltzmann2DCheckerboard<E, D, d, R, momentOrder> ProblemType;
  static const size_t dimRange     = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Tests::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  Boltzmann2DCheckerboardTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(divide_t_end_by, Stuff::Grid::Providers::Cube<G>::create(ProblemType::default_grid_config())->grid_ptr(),
               num_refs)
    , problem_(*(ProblemType::create(ProblemType::default_config())))
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual bool provides_exact_solution() const override final
  {
    return false;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+======================================================================================================+\n"
        << "|+====================================================================================================+|\n"
        << "||  Testcase: Boltzmann 2D Checkerboard                                                               ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [0, 7] x [0, 7]                                                                          ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end())
               + "]                                                                                  ||\n"
        << "||  flux = see http://dx.doi.org/10.1016/j.jcp.2005.04.011 Section 4.1                                ||\n"
        << "||  rhs = see http://dx.doi.org/10.1016/j.jcp.2005.04.011 Section 4.1                                 ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
}; // class Boltzmann2DCheckerboardTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH
