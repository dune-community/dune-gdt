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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/la/container.hh>

#include "../default.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** Testcase for the \f$P_n\f$ moment approximation of the Fokker-Planck equation in one dimension
 * \f[
 * \partial_t \psi(t,x,v) + v * \partial_x \psi(t,x,v) + \sigma_a(x)*\psi(t,x,v) = 0.5*T(x)*\Delta_v \psi(t,x,v) + Q(x),
 * \f]
 * where \f$\psi: [0,T] \times [x_l, x_r] \times [-1, 1] \to \mathbb{R}\f$,
 * \f$Delta_v \psi = \partial_v( (1-v^2)\partial_v \psi)\f$ is the Laplace-Beltrami operator and
 * \f$\sigma_a, T, Q: [x_l, x_r] \to \mathbb{R}\f$ are the absorption coefficient, the transport coefficient and a
 * source, respectively.
 * The \f$P_n\f$ model approximates the solution of the Fokker-Planck equation by an ansatz
 * \f[
 * \psi(t,x,v) = \sum \limits_{l=0}^n u_i(t,x)\phi_i(v)
 * \f]
 * where the \f$\phi_i, \, i=0,\ldots,n$ are suitable basis functions of (a subset of) the function space on
 * \f$[-1, 1]\f$ that \f$\psi\f$ lives in. n is called the moment order. Usually, the \f$\phi_i\f$ are chosen as the
 * Legendre polynomials up to order n.
 * Once suitable basis functions are found, a Galerkin semidiscretization in v is done, so the \f$\phi_i\f$ are also
 * taken as basis for the test space. This results in an equation of the form
 * \f[
 * M \partial_t u + D \partial_x u = q - (\sigma_a*M + 0.5*T*S) u,
 * \f]
 *  where \f$u = (u_1, \ldots, u_n)^T\f$, \f$M, D, S \in \mathbb{R}^{n\times n}\f$ with
 * \f$ M_{ji} = (\phi_i, \phi_j)_v\f$, \f$D_{ji} = (v*\phi_i, \phi_j)_v\f$,
 * \f$S_{ji} = ((1-v^2)\partial_v \phi_i, \partial_v \phi_j)_v\f$ and
 * \f$q_i(x) = (Q(x), \phi_i(v))_v = Q(x) (1, \phi_i(v))_v\f$. Here, \f$(a,b)_v = \int \limits_{-1}^1 a(v)b(v)dv\f$
 * is the \f$L^2\f$ inner product with respect to \f$v\f$.
 * In the following, we rescale \f$u\f$ s.t. \f$(\psi(t,x,v),\phi_i(v))_v = u_i(t,x)\f$ if the ansatz holds. Provided
 * the \f$ \phi_i \f$ are an orthogonal basis, \f$M\f$ is invertible and the rescaling corresponds to a multiplication
 * of \f$u\f$ by \f$M^{-1}\f$ from the left, giving the equation
 * \f[
 * \partial_t u + D M^{-1} \partial_x u = q - (\sigma_a*I_{n\times n} + 0.5*T*S M^{-1}) u.
 * \f]
 * For details on the parameters of the test cases implemented here (OneBeam, TwoBeams, TwoPulses, RectangularIC,
 * SourceBeam) see Schneider, Alldredge, Frank, Klar, "Higher Order Mixed-Moment Approximations for the
 * Fokker-Planck Equation in One Space Dimension", SIAM J. Appl. Math., 74(4), 1087–1114
 */
template <class E, class D, size_t d, class R, size_t momentOrder>
class TwoBeamsPnLegendre : public Default<E, D, d, R, momentOrder + 1>
{
  typedef TwoBeamsPnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef Default<E, D, d, R, momentOrder + 1> BaseType;

public:
  static const bool linear = true;
  static const int precision = 15; // precision for to_string
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::DummyEntityType;
  typedef typename Dune::XT::Functions::AffineFunction<DummyEntityType, R, dimRange, R, dimRange, dimDomain>
      FluxAffineFunctionType;
  typedef typename Dune::GDT::GlobalFunctionBasedAnalyticalFlux<FluxAffineFunctionType, E, D, d, R, dimRange, 1>
      DefaultFluxType;
  typedef typename DefaultFluxType::RangeType RangeType;
  typedef typename DefaultFluxType::FluxRangeType FluxRangeType;
  typedef typename FluxAffineFunctionType::FieldMatrixType MatrixType;
  using typename BaseType::DefaultInitialValueType;
  typedef typename XT::Functions::AffineFunction<DummyEntityType, R, dimRange, R, dimRange, 1> RHSAffineFunctionType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<RHSAffineFunctionType, E, D, d, R, dimRange, 1>
      RHSCheckerboardFunctionType;
  typedef typename Dune::GDT::CheckerboardBasedRhsEvaluation<RHSCheckerboardFunctionType, E, D, d, R, dimRange, 1>
      DefaultRHSType;
  typedef typename DefaultRHSType::DomainType DomainType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".twobeams_Pn_legendre";
  }

  std::string type() const override
  {
    return BaseType::type() + ".twobeams_Pn_legendre";
  }

  static std::string short_id()
  {
    return "2BeamsPnLegendre";
  }

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[1.0]";
    grid_config["num_elements"] = "[100]";
    grid_config["overlap_size"] = "[1 1 1 1]";
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
    const std::shared_ptr<const DefaultRHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const DefaultInitialValueType> initial_values(
        DefaultInitialValueType::create(config.sub("initial_values")));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const DefaultBoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const double psi_vac = 1e-4)
  {
    ConfigType config;
    config.add(default_grid_config(), "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    flux_config["A"] = create_flux_matrix();
    flux_config["b"] = Dune::XT::Common::to_string(RangeType(0));
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"] = "[0.0]";
    rhs_config["upper_right"] = "[1.0]";
    rhs_config["num_elements"] = "[1]";
    create_rhs_values(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[1.0]";
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = create_initial_values(psi_vac);
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"] = BoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = create_boundary_values();
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values");
    return config;
  } // ... default_config(...)

  TwoBeamsPnLegendre(const std::shared_ptr<const FluxType> flux_in,
                     const std::shared_ptr<const RHSType> rhs_in,
                     const std::shared_ptr<const InitialValueType> initial_values_in,
                     const ConfigType& grid_config_in,
                     const ConfigType& boundary_info_in,
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
    return 4.0;
  }

  virtual bool has_non_zero_rhs() const override
  {
    return true;
  }

protected:
  // n-th component of RHS is (T/2 n(n+1) - sigma_a) u_n + 2Q delta(n).
  // here, T = 0, Q = 0, sigma_a = 4
  static void create_rhs_values(ConfigType& rhs_config)
  {
    Dune::FieldMatrix<R, dimRange, dimRange> A(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      A[rr][rr] = -4;
    rhs_config["A.0"] = XT::Common::to_string(A, precision);
    rhs_config["b.0"] = Dune::XT::Common::to_string(FluxRangeType(0));
  } // ... create_rhs_values(...)

  // flux matrix F[rr][cc] = rr/(2*rr + 1)       if cc == rr - 1
  //                       = (rr + 1)/(2*rr + 1) if cc == rr + 1
  //                       = 0                   else
  static std::string create_flux_matrix()
  {
    Dune::FieldMatrix<R, dimRange, dimRange> A(0);
    for (size_t rr = 0; rr < dimRange; ++rr) {
      for (size_t cc = 0; cc < dimRange; ++cc) {
        if (cc == rr - 1)
          A[rr][cc] = rr / (2.0 * rr + 1.0);
        else if (cc == rr + 1)
          A[rr][cc] = (rr + 1.) / (2.0 * rr + 1.0);
      }
    }
    return XT::Common::to_string(A, precision);
  } // ... create_flux_matrix()

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is 2 psi_vac if n == 0 and 0 else.
  static std::string create_initial_values(const R psi_vac = 1e-4)
  {
    Dune::FieldVector<R, dimRange> initial_vals(0);
    initial_vals[0] = 2 * psi_vac;
    return XT::Common::to_string(initial_vals, precision);
  } // ... create_initial_values()

  // boundary value of kinetic equation is 100*delta(v-1) at x = 0 and 100*delta(v+1) at x = 1,
  // so k-th component of boundary value has to be 50*\phi_k(1) at x = 0 and 50*\phi_k(-1) at x = 1.
  // Model with function(x) = 50*((\phi_k(-1) - \phi_k(1))*x + \phi_k(1)).
  // For Legendre polynomials, this is [50 50 50 ...] at x = 0 and [50 -50 50 -50 ... ] at x = 1,
  // so the function used is function(x) = 50*((-1)^n - 1)*x + 1)
  static std::string create_boundary_values()
  {
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += "50*(" + Dune::XT::Common::to_string(((1.0 - 2.0 * (rr % 2)) - 1.0), precision) + "*x[0]+1)";
    }
    str += "]";
    return str;
  } // ... create_boundary_values()

}; // class TwoBeamsPnLegendre<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_TWOBEAMS_HH
