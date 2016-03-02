// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH
#define DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/**
 * Testcase for the Boltzmann equation in two dimensions,
 * see Brunner, Holloway, "Two-dimensional time dependent Riemann solvers for neutron transport", Journal of
 * Computational Physics, Volume 210, Issue 1, 2005
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
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DummyEntityType;
  typedef typename Dune::Stuff::Functions::Affine<DummyEntityType, RangeFieldImp, dimRange, RangeFieldImp, dimRange,
                                                  dimDomain> FluxAffineFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxAffineFunctionType, EntityImp, DomainFieldImp,
                                                                dimDomain, RangeFieldImp, dimRange, 1> DefaultFluxType;
  typedef typename DefaultFluxType::RangeType FluxRangeType;
  typedef typename DefaultFluxType::MatrixType MatrixType;
  using typename BaseType::DefaultInitialValueType;
  typedef typename DS::Functions::Affine<DummyEntityType, RangeFieldImp, dimRange, RangeFieldImp, dimRange, 1>
      RHSAffineFunctionType;
  typedef typename DS::Functions::FunctionCheckerboard<RHSAffineFunctionType, EntityImp, DomainFieldImp, dimDomain,
                                                       RangeFieldImp, dimRange, 1> RHSCheckerboardFunctionType;
  typedef typename Dune::GDT::CheckerboardBasedRHS<RHSCheckerboardFunctionType, EntityImp, DomainFieldImp, dimDomain,
                                                   RangeFieldImp, dimRange, 1> DefaultRHSType;
  typedef typename DefaultRHSType::RangeRangeType RangeType;
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
    static const int precision = 15; // precision for toString

    // (\Sigma_s \delta_{l0}\delta{m0} - \Sigma_t) * \psi_l^m
    static void create_source_values(ConfigType& source_config)
    {
      source_config["lower_left"]   = "[0.0 0.0]";
      source_config["upper_right"]  = "[1.0 1.0]";
      source_config["num_elements"] = "[1 1]";
      const RangeFieldImp Sigma_s   = 0;
      const RangeFieldImp Sigma_t   = 0;
      MatrixType S;
      S[pos(0, 0)][pos(0, 0)] = Sigma_s - Sigma_t;
      for (size_t l = 1; l <= momentOrder; ++l)
        for (size_t m = 0; m <= l; ++m)
          S[pos(l, m)][pos(l, m)] = -1.0 * Sigma_t;
      source_config["A.0"]      = DSC::toString(S, precision);
      source_config["b.0"]      = DSC::toString(RangeType(0));
      source_config["sparse.0"] = "true";
    } // ... create_source_values(...)

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
            Z[row][pos(l + 1, m)] = B(l + 1, m);
          }
        }
      }
      flux_config["A.0"]      = DSC::toString(X, precision);
      flux_config["sparse.0"] = "true";
      flux_config["A.1"]      = DSC::toString(Z, precision);
      flux_config["sparse.1"] = "true";
      flux_config["b"]        = DSC::toString(typename DefaultFluxType::RangeType(0));
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
  // Domain should be [-0.5, 0.5]^2, change it to [0, 1]^2 to make YaspGrid happy
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

  static std::unique_ptr<ThisType> create(const ConfigType cfg = default_config(),
                                          const std::string sub_name = static_id())
  {
    const ConfigType config = cfg.has_sub(sub_name) ? cfg.sub(sub_name) : cfg;
    const std::shared_ptr<const DefaultFluxType> flux(DefaultFluxType::create(config.sub("flux")));
    RangeType alpha;
    alpha[0] = std::log(2);
    //    const std::shared_ptr< const DefaultFluxType > flux
    //        = std::make_shared< const DefaultFluxType >(GetData::velocity_grid_view(),
    //                                                    GetData::basefunctions(),
    //                                                    alpha,
    //                                                    0.5,
    //                                                    10e-8,
    //                                                    0.01,
    //                                                    0.001);
    const std::shared_ptr<const DefaultRHSType> source(DefaultRHSType::create(config.sub("source")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config   = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return Stuff::Common::make_unique<ThisType>(
        flux, source, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const std::string sub_name = "")
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    GetData::create_flux_matrices(flux_config);
    config.add(flux_config, "flux");
    ConfigType source_config;
    ;
    GetData::create_source_values(source_config);
    config.add(source_config, "source");
    ConfigType initial_value_config      = DefaultInitialValueType::default_config();
    initial_value_config["lower_left"]   = "[0.0 0.0]";
    initial_value_config["upper_right"]  = "[1.0 1.0]";
    initial_value_config["num_elements"] = "[1 1]";
    initial_value_config["variable"]     = "x";
    initial_value_config["values.0"]     = GetData::create_initial_values();
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config    = BoundaryValueType::default_config();
    boundary_value_config["type"]       = BoundaryValueType::static_id();
    boundary_value_config["variable"]   = "x";
    boundary_value_config["expression"] = GetData::create_boundary_values();
    boundary_value_config["order"] = "0";
    config.add(boundary_value_config, "boundary_values");
    if (sub_name.empty())
      return config;
    else {
      ConfigType tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  Boltzmann2DLineSource(const std::shared_ptr<const FluxType> flux_in, const std::shared_ptr<const RHSType> source_in,
                        const std::shared_ptr<const InitialValueType> initial_values_in,
                        const ConfigType& grid_config_in, const ConfigType& boundary_info_in,
                        const std::shared_ptr<const BoundaryValueType> boundary_values_in)
    : BaseType(flux_in, source_in, initial_values_in, grid_config_in, boundary_info_in, boundary_values_in)
  {
  }
}; // ... Boltzmann2DLineSource ...


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_HDD_HYPERBOLIC_PROBLEMS_2DBOLTZMANN_HH
