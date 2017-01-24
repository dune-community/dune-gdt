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

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/xt/common/string.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

double binomial_coefficient(const double n, const size_t k)
{
  double ret(1);
  for (size_t ii = 1; ii <= k; ++ii)
    ret *= (n + 1 - ii) / ii;
  return ret;
}

template <class FieldType>
FieldType evaluate_legendre_polynomial(const FieldType& v, const size_t n)
{
  FieldType ret(0);
  for (size_t k = 0; k <= n; ++k)
    ret += std::pow(-v, k) * binomial_coefficient(n, k) * binomial_coefficient((n + k - 1) / 2., n);
  ret *= std::pow(-1, n) * (1 << n); // 2^n
  return ret;
}

template <class FieldType, class RangeType>
FieldType evaluate_hat_function(const FieldType& v, const size_t n, const RangeType& v_points)
{
  if (n == 0 && XT::Common::FloatCmp::ge(v, v_points[0]) && XT::Common::FloatCmp::le(v, v_points[1]))
    return (v - v_points[1]) / (v_points[0] - v_points[1]);
  else if (n == v_points.size() && XT::Common::FloatCmp::ge(v, v_points[v_points.size() - 2])
           && XT::Common::FloatCmp::le(v, v_points[v_points.size() - 1]))
    return (v - v_points[v_points.size() - 2]) / (v_points[v_points.size() - 1] - v_points[v_points.size() - 2]);
  else if (XT::Common::FloatCmp::ge(v, v_points[n - 1]) && XT::Common::FloatCmp::le(v, v_points[n + 1]))
    if (XT::Common::FloatCmp::le(v, v_points[n]))
      return (v - v_points[n - 1]) / (v_points[n] - v_points[n - 1]);
    else
      return (v - v_points[n + 1]) / (v_points[n] - v_points[n + 1]);
  else
    return 0;
}

template <class FieldType, class PointsVectorType>
FieldType evaluate_first_order_dg(const FieldType& v, const size_t n, const PointsVectorType& v_points)
{
  if (XT::Common::FloatCmp::ge(v, v_points[n / 2]) && XT::Common::FloatCmp::le(v, v_points[n / 2 + 1]))
    if (n % 2)
      return v;
    else
      return 1;
  else
    return 0;
}


/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
class SourceBeamPnLegendre : public TwoBeamsPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder>
{
  typedef SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> ThisType;
  typedef TwoBeamsPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".sourcebeam_Pn_legendre";
  }

  std::string type() const override
  {
    return BaseType::type() + ".sourcebeam_Pn_legendre";
  }

  static std::string short_id()
  {
    return "SourceBeamPnLegendre";
  }

public:
  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0.0]";
    grid_config["upper_right"] = "[3.0]";
    grid_config["num_elements"] = "[3000]";
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

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const bool BGK_type_collision = true,
                                   const RangeFieldImp psi_vac = 5e-9)
  {
    ConfigType config = BaseType::default_config(psi_vac);
    config.add(grid_config, "grid", true);
    config.add(default_boundary_info_config(), "boundary_info", true);
    ConfigType rhs_config;
    rhs_config["lower_left"] = "[0.0]";
    rhs_config["upper_right"] = "[3.0]";
    rhs_config["num_elements"] = "[6]";
    BGK_type_collision ? create_rhs_values_BGK(rhs_config) : create_rhs_values_laplace_beltrami(rhs_config);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs", true);
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[3.0]";
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = create_initial_values(psi_vac);
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values", true);
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] =
        BGK_type_collision ? create_realizable_boundary_values(psi_vac) : create_boundary_values();
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values", true);
    return config;
  } // ... default_config(...)

  SourceBeamPnLegendre(const std::shared_ptr<const FluxType> flux_in,
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
  using BaseType::create_initial_values;

  // n-th component of RHS is -(T/2 n(n+1) + sigma_a) u_n + 2Q delta(n)
  // For this test case,
  // sigma_a = 1 if x <= 2, 0 else
  // T = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  // Thus,
  // rhs[n] = -u[n]                                     if x <= 1
  //        = -(1 + n(n+1))u[n] + 2*delta(n)            if 1 < x <= 1.5
  //        = -(1 + n(n+1))u[n]                         if 1.5 < x <= 2
  //        = -5*n(n+1)*u[n]                            else (2 < x <= 3)
  static void create_rhs_values_laplace_beltrami(ConfigType& rhs_config)
  {
    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldImp, dimRange> b(0);
      if (ii == 2) // 1 < x <= 1.5
        b[0] = 2;
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (ii == 0 || ii == 1) // x <= 1
          A[rr][rr] = -1;
        else if (ii == 2 || ii == 3) // 1 <= x <= 2
          A[rr][rr] = -1. - rr * (rr + 1);
        else // 2 <= x <= 3
          A[rr][rr] = -5. * rr * (rr + 1);
      } // rr
      rhs_config["A." + XT::Common::to_string(ii)] = XT::Common::to_string(A);
      rhs_config["b." + XT::Common::to_string(ii)] = XT::Common::to_string(b);
    } // ii
  } // ... create_rhs_values_laplace_beltrami()

  // n-th component of RHS is -sigma_t u_n + (sigma_s u_0 + 2Q) delta(n)
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  // Thus,
  // rhs[n] = -u[n]                                     if x <= 1
  //        = (2 delta(n) - 3)u[n] + 2*delta(n)         if 1 < x <= 1.5
  //        = (2 delta(n) - 3)u[n]                      if 1.5 < x <= 2
  //        = (10 delta(n) - 10)*u[n]                   else (2 < x <= 3)
  static void create_rhs_values_BGK(ConfigType& rhs_config)
  {
    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldImp, dimRange> b(0);
      if (ii == 2) // 1 < x <= 1.5
        b[0] = 2;
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (ii == 0 || ii == 1) // x <= 1
          A[rr][rr] = -1;
        else if (ii == 2 || ii == 3) // 1 <= x <= 2
          A[rr][rr] = 2 * (rr == 0) - 3;
        else // 2 <= x <= 3
          A[rr][rr] = 10 * (rr == 0) - 10;
      } // rr
      rhs_config["A." + XT::Common::to_string(ii)] = XT::Common::to_string(A);
      rhs_config["b." + XT::Common::to_string(ii)] = XT::Common::to_string(b);
    } // ii
  } // ... create_rhs_values_BGK()

  // Boundary value of kinetic equation is delta(v-1) at x = 0 and psi_vac at x = 3,
  // so n-th component of boundary value has to be 0.5*\phi_n(1) at x = 0 and 0.5*psi_vac*delta(n) at x = 3.
  // For Legendre polynomials, this is [0.5 0.5 0.5 ...] at x = 0 and [2*psi_vac 0 0 ... ] at x = 3.
  // Model with linear interpolating function.
  static std::string create_boundary_values(const RangeFieldImp psi_vac = 1e-4)
  {
    std::string str = "[0.5-" + XT::Common::to_string(0.5 - 2 * psi_vac, precision) + "*x[0]/3.0";
    for (size_t cc = 1; cc < dimRange; ++cc)
      str += " 0.5-0.5*x[0]/3.0";
    str += "]";
    return str;
  } // ... create_boundary_values()

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // For Legendre polynomials, this is approx. [1 0.98 0.98 ...] at x = 0 and [2*\psi_{vac} 0 0 ... ]
  // at x = 3.
  // Simulate with linear interpolating function.
  static std::string create_realizable_boundary_values(const RangeFieldImp psi_vac = 5e-9)
  {
    std::string str = "[" + XT::Common::to_string(left_boundary_vals[0], precision) + "-("
                      + XT::Common::to_string(left_boundary_vals[0] - 2 * psi_vac, precision) + ")*x[0]/3.0";
    for (size_t cc = 1; cc < dimRange; ++cc)
      str += " " + XT::Common::to_string(left_boundary_vals[cc], precision) + "-("
             + XT::Common::to_string(left_boundary_vals[cc], precision) + ")*x[0]/3.0";
    str += "]";
    return str;
  } // ... create_realizable_boundary_values()

private:
  static std::vector<RangeFieldImp> left_boundary_vals;
};

template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t momentOrder>
std::vector<RangeFieldImp>
    SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, momentOrder>::left_boundary_vals = {
        1.0,
        0.99821587588384722888546103381424,
        0.99465512765154168665638310109402,
        0.98933271069998046949348833831956,
        0.98227094694487696211966529188208,
        0.973499392321142130859194712695,
        0.96305466163489000577549441482268,
        0.95098021205162953191200986872772,
        0.93732608680563406367139508178687,
        0.9221486210023942981265256212649,
        0.9055101116570464980458926840956,
        0.88747845436446410855581951816166,
        0.86812674922921981743467893406176,
        0.84753287889397251365597836272049,
        0.82577906169132591028268357189113,
        0.80295138310538515338857079818232,
        0.77913930886388561354086750071819,
        0.75443518308891443080916196222167,
        0.72893371501317733202250980561682,
        0.70273145781902986912513049892392,
        0.67592628317890714138498882739631,
        0.6486168550684271196351824164262,
        0.62090210638764909412176999841986,
        0.5928807218623355455816081929013,
        0.56465063060643566894968702250824,
        0.53630851161046188339453821024199,
        0.50794931527927140373971781871463,
        0.47966580397850727702579631268181,
        0.4515481143633037274085705597037,
        0.42368334405769961332718708268865,
        0.39615516503056535291425143183468,
        0.3690434657758989934336017941775};


/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t num_points>
class SourceBeamPnHatFunctions
    : public SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points - 1>
{
  typedef SourceBeamPnHatFunctions<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points> ThisType;
  typedef SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points - 1> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".sourcebeam_Pn_hatfunctions";
  }

  std::string type() const override
  {
    return BaseType::type() + ".sourcebeam_Pn_hatfunctions";
  }

  static std::string short_id()
  {
    return "SourceBeamPnHatFunctions";
  }

  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;

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

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const RangeType v_points = create_equidistant_points(),
                                   const RangeFieldImp psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    flux_config["A"] = create_flux_matrix(v_points);
    flux_config["b"] = Dune::XT::Common::to_string(RangeType(0));
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"] = "[0.0]";
    rhs_config["upper_right"] = "[3.0]";
    rhs_config["num_elements"] = "[6]";
    create_rhs_values(rhs_config, v_points);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[3.0]";
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = create_initial_values(v_points, psi_vac);
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = create_boundary_values(v_points, psi_vac);
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values");
    return config;
  } // ... default_config(...)

  SourceBeamPnHatFunctions(const std::shared_ptr<const FluxType> flux_in,
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

  static RangeType create_equidistant_points()
  {
    RangeType ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = -1. + 2. * ii / (dimRange - 1.);
    return ret;
  }

protected:
  // returns <b>, where b is the basis functions vector
  static RangeType hatfunctions_integrated(const RangeType& v_points)
  {
    RangeType ret(0);
    ret[0] = v_points[1] - v_points[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      ret[ii] = v_points[ii + 1] - v_points[ii - 1];
    ret[dimRange - 1] = v_points[dimRange - 1] - v_points[dimRange - 2];
    ret *= 0.5;
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  static MatrixType mass_matrix(const RangeType& v_points)
  {
    MatrixType ret(0);
    ret[0][0] = (v_points[1] - v_points[0]) / 3.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (v_points[rr + 1] - v_points[rr - 1]) / 3.;
      if (rr > 0)
        ret[rr][rr - 1] = (v_points[rr] - v_points[rr - 1]) / 6.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] = (v_points[rr + 1] - v_points[rr]) / 6.;
    }
    ret[dimRange - 1][dimRange - 1] = (v_points[dimRange - 1] - v_points[dimRange - 2]) / 3.;
    return ret;
  }

  // returns matrix with entries <v h_i h_j>
  static MatrixType mass_matrix_with_v(const RangeType& v_points)
  {
    MatrixType ret(0);
    ret[0][0] = (v_points[1] * v_points[1] + 2 * v_points[1] * v_points[0] - 3 * v_points[0] * v_points[0]) / 12.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (v_points[rr + 1] * v_points[rr + 1] + 2 * v_points[rr + 1] * v_points[rr]
                       - 2 * v_points[rr] * v_points[rr - 1]
                       - v_points[rr - 1] * v_points[rr - 1])
                      / 12.;
      if (rr > 0)
        ret[rr][rr - 1] = (v_points[rr] * v_points[rr] - v_points[rr - 1] * v_points[rr - 1]) / 12.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] = (v_points[rr + 1] * v_points[rr + 1] - v_points[rr] * v_points[rr]) / 12.;
    }
    ret[dimRange - 1][dimRange - 1] =
        (3 * v_points[dimRange - 1] * v_points[dimRange - 1] - 2 * v_points[dimRange - 1] * v_points[dimRange - 2]
         - v_points[dimRange - 2] * v_points[dimRange - 2])
        / 12.;
    return ret;
  }

  // see https://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
  static MatrixType tridiagonal_matrix_inverse(const MatrixType& matrix)
  {
#ifndef NDEBUG
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        if ((cc > rr + 1 || cc + 1 < rr) && XT::Common::FloatCmp::ne(matrix[rr][cc], 0.))
          DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Matrix has to be tridiagonal!");
#endif // NDEBUG
    MatrixType ret(0);
    Dune::FieldVector<RangeFieldImp, dimRange + 1> a(0), b(0), c(0), theta(0);
    Dune::FieldVector<RangeFieldImp, dimRange + 2> phi(0);
    for (size_t ii = 1; ii < dimRange + 1; ++ii) {
      a[ii] = matrix[ii - 1][ii - 1];
      if (ii < dimRange) {
        b[ii] = matrix[ii - 1][ii];
        c[ii] = matrix[ii][ii - 1];
      }
    }
    theta[0] = 1;
    theta[1] = a[1];
    for (size_t ii = 2; ii < dimRange + 1; ++ii)
      theta[ii] = a[ii] * theta[ii - 1] - b[ii - 1] * c[ii - 1] * theta[ii - 2];
    phi[dimRange + 1] = 1;
    phi[dimRange] = a[dimRange];
    for (size_t ii = dimRange - 1; ii > 0; --ii)
      phi[ii] = a[ii] * phi[ii + 1] - b[ii] * c[ii] * phi[ii + 2];
    for (size_t ii = 1; ii < dimRange + 1; ++ii) {
      for (size_t jj = 1; jj < dimRange + 1; ++jj) {
        if (ii == jj)
          ret[ii - 1][jj - 1] = theta[ii - 1] * phi[jj + 1] / theta[dimRange];
        else if (ii < jj) {
          ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[ii - 1] * phi[jj + 1] / theta[dimRange];
          for (size_t kk = ii; kk < jj; ++kk)
            ret[ii - 1][jj - 1] *= b[kk];
        } else if (ii > jj) {
          ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[jj - 1] * phi[ii + 1] / theta[dimRange];
          for (size_t kk = jj; kk < ii; ++kk)
            ret[ii - 1][jj - 1] *= c[kk];
        }
      } // jj
    } // ii
    return ret;
  } // ... tridiagonal_matrix_inverse(...)

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  static std::string create_flux_matrix(const RangeType& v_points)
  {
    MatrixType A = mass_matrix_with_v(v_points);
    const MatrixType M = mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    A.rightmultiply(M_inv);
    return XT::Common::to_string(A, precision);
  } // ... create_flux_matrix

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the moment vector is psi_vac * <b>.
  static std::string create_initial_values(const RangeType& v_points = create_equidistant_points(),
                                           const RangeFieldImp psi_vac = 5e-9)
  {
    RangeType initial_vals = hatfunctions_integrated(v_points);
    initial_vals *= psi_vac;
    return XT::Common::to_string(initial_vals, precision);
  } // ... create_initial_values()

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  static void create_rhs_values(ConfigType& rhs_config, const RangeType& v_points)
  {
    const RangeType integrated_basis = hatfunctions_integrated(v_points);
    const MatrixType M = mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);

    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldImp, dimRange> b(0);
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (ii == 2) // 1 < x <= 1.5
          b[rr] = integrated_basis[rr];
        if (ii == 0 || ii == 1) // x <= 1
          A[rr][rr] = -1;
        for (size_t cc = 0; cc < dimRange; ++cc) {
          if (ii == 2 || ii == 3) // 1 <= x <= 2
            A[rr][cc] = integrated_basis[rr] * c[cc] - 3. * (rr == cc);
          else if (ii == 4 || ii == 5) // 2 <= x <= 3
            A[rr][cc] = 5 * integrated_basis[rr] * c[cc] - 10. * (rr == cc);
        } // cc
      } // rr
      rhs_config["A." + XT::Common::to_string(ii)] = XT::Common::to_string(A);
      rhs_config["b." + XT::Common::to_string(ii)] = XT::Common::to_string(b);
    } // ii
  } // ... create_rhs_values()

  // returns the numerator g of the left boundary value (see create_boundary_values)
  static RangeFieldImp numerator(const RangeFieldImp v)
  {
    return std::exp(-1e5 * (v - 1) * (v - 1));
  }

  // returns the denominator <g> of the left boundary value (see create_boundary_values)
  static RangeFieldImp denominator()
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10) * std::erf(200 * std::sqrt(10));
  }

  // calculates integral from v_l to v_u of numerator g
  static RangeFieldImp integral_1(RangeFieldImp v_l, RangeFieldImp v_u)
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10)
           * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
  }

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // Simulate with linear interpolating function.
  static std::string create_boundary_values(const RangeType& v_points = create_equidistant_points(),
                                            const RangeFieldImp psi_vac = 5e-9)
  {
    const RangeType integrated_basis = hatfunctions_integrated(v_points);
    RangeType left_boundary_values(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      const auto vnm = v_points[nn - 1];
      const auto vn = v_points[nn];
      const auto vnp = v_points[nn + 1];
      if (nn < dimRange - 1)
        left_boundary_values[nn] += 1. / ((vn - vnp) * denominator())
                                    * ((1 - vnp) * integral_1(vn, vnp) - 1. / 2e5 * (numerator(vnp) - numerator(vn)));
      if (nn > 0)
        left_boundary_values[nn] += 1. / ((vn - vnm) * denominator())
                                    * ((1 - vnm) * integral_1(vnm, vn) - 1. / 2e5 * (numerator(vn) - numerator(vnm)));
    }

    std::string str = "[";
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0)
        str += " ";
      str += XT::Common::to_string(left_boundary_values[nn], precision) + "-("
             + XT::Common::to_string(left_boundary_values[nn] - psi_vac * integrated_basis[nn], precision)
             + ")*x[0]/3.0";
    }
    str += "]";
    return str;
  } // ... create_boundary_values()
}; // class SourceBeamPnHatFunctions


/** \see class TwoBeams in twobeams.hh */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t num_points>
class SourceBeamPnFirstOrderDG
    : public SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 2 * num_points - 3>
{
  typedef SourceBeamPnFirstOrderDG<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, num_points> ThisType;
  typedef SourceBeamPnLegendre<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 2 * num_points - 3> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::FluxType;
  using typename BaseType::RHSType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;

  typedef typename Dune::FieldVector<DomainFieldImp, num_points> PointsVectorType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".sourcebeam_Pn_firstorderdg";
  }

  std::string type() const override
  {
    return BaseType::type() + ".sourcebeam_Pn_firstorderdg";
  }

  static std::string short_id()
  {
    return "SourceBeamPnFirstOrderDG";
  }

  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;

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

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const PointsVectorType& v_points = create_equidistant_points(),
                                   const RangeFieldImp psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    flux_config["A"] = create_flux_matrix(v_points);
    flux_config["b"] = Dune::XT::Common::to_string(RangeType(0));
    config.add(flux_config, "flux");
    ConfigType rhs_config;
    rhs_config["lower_left"] = "[0.0]";
    rhs_config["upper_right"] = "[3.0]";
    rhs_config["num_elements"] = "[6]";
    create_rhs_values(rhs_config, v_points);
    rhs_config["name"] = static_id();
    config.add(rhs_config, "rhs");
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = "[0.0]";
    initial_value_config["upper_right"] = "[3.0]";
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["values.0"] = create_initial_values(v_points, psi_vac);
    initial_value_config["name"] = static_id();
    config.add(initial_value_config, "initial_values");
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["expression"] = create_boundary_values(v_points, psi_vac);
    boundary_value_config["order"] = "10";
    config.add(boundary_value_config, "boundary_values");
    return config;
  } // ... default_config(...)

  SourceBeamPnFirstOrderDG(const std::shared_ptr<const FluxType> flux_in,
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

  static PointsVectorType create_equidistant_points()
  {
    PointsVectorType ret;
    for (size_t ii = 0; ii < num_points; ++ii)
      ret[ii] = -1. + 2. * ii / (num_points - 1);
    return ret;
  }

protected:
  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& v_points)
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < num_points - 1; ++ii) {
      ret[2 * ii] = v_points[ii + 1] - v_points[ii];
      ret[2 * ii + 1] = (std::pow(v_points[ii + 1], 2) - std::pow(v_points[ii], 2)) / 2.;
    }
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  static MatrixType mass_matrix(const PointsVectorType& v_points)
  {
    MatrixType ret(0);
    for (size_t ii = 0; ii < num_points - 1; ++ii) {
      ret[2 * ii][2 * ii] = v_points[ii + 1] - v_points[ii];
      ret[2 * ii + 1][2 * ii + 1] = (std::pow(v_points[ii + 1], 3) - std::pow(v_points[ii], 3)) / 3.;
      ret[2 * ii][2 * ii + 1] = (std::pow(v_points[ii + 1], 2) - std::pow(v_points[ii], 2)) / 2.;
      ret[2 * ii + 1][2 * ii] = ret[2 * ii][2 * ii + 1];
    }
    return ret;
  }

  // returns matrix with entries <v h_i h_j>
  static MatrixType mass_matrix_with_v(const PointsVectorType& v_points)
  {
    MatrixType ret(0);
    for (size_t ii = 0; ii < num_points - 1; ++ii) {
      ret[2 * ii][2 * ii] = (std::pow(v_points[ii + 1], 2) - std::pow(v_points[ii], 2)) / 2.;
      ret[2 * ii + 1][2 * ii + 1] = (std::pow(v_points[ii + 1], 4) - std::pow(v_points[ii], 4)) / 4.;
      ret[2 * ii][2 * ii + 1] = (std::pow(v_points[ii + 1], 3) - std::pow(v_points[ii], 3)) / 3.;
      ret[2 * ii + 1][2 * ii] = ret[2 * ii][2 * ii + 1];
    }
    return ret;
  }

  // see https://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
  static MatrixType tridiagonal_matrix_inverse(const MatrixType& matrix)
  {
#ifndef NDEBUG
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        if ((cc > rr + 1 || cc + 1 < rr) && XT::Common::FloatCmp::ne(matrix[rr][cc], 0.))
          DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Matrix has to be tridiagonal!");
#endif // NDEBUG
    MatrixType ret(0);
    Dune::FieldVector<RangeFieldImp, dimRange + 1> a(0), b(0), c(0), theta(0);
    Dune::FieldVector<RangeFieldImp, dimRange + 2> phi(0);
    for (size_t ii = 1; ii < dimRange + 1; ++ii) {
      a[ii] = matrix[ii - 1][ii - 1];
      if (ii < dimRange) {
        b[ii] = matrix[ii - 1][ii];
        c[ii] = matrix[ii][ii - 1];
      }
    }
    theta[0] = 1;
    theta[1] = a[1];
    for (size_t ii = 2; ii < dimRange + 1; ++ii)
      theta[ii] = a[ii] * theta[ii - 1] - b[ii - 1] * c[ii - 1] * theta[ii - 2];
    phi[dimRange + 1] = 1;
    phi[dimRange] = a[dimRange];
    for (size_t ii = dimRange - 1; ii > 0; --ii)
      phi[ii] = a[ii] * phi[ii + 1] - b[ii] * c[ii] * phi[ii + 2];
    for (size_t ii = 1; ii < dimRange + 1; ++ii) {
      for (size_t jj = 1; jj < dimRange + 1; ++jj) {
        if (ii == jj)
          ret[ii - 1][jj - 1] = theta[ii - 1] * phi[jj + 1] / theta[dimRange];
        else if (ii < jj) {
          ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[ii - 1] * phi[jj + 1] / theta[dimRange];
          for (size_t kk = ii; kk < jj; ++kk)
            ret[ii - 1][jj - 1] *= b[kk];
        } else if (ii > jj) {
          ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[jj - 1] * phi[ii + 1] / theta[dimRange];
          for (size_t kk = jj; kk < ii; ++kk)
            ret[ii - 1][jj - 1] *= c[kk];
        }
      } // jj
    } // ii
    return ret;
  } // ... tridiagonal_matrix_inverse(...)

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  static std::string create_flux_matrix(const PointsVectorType& v_points)
  {
    MatrixType A = mass_matrix_with_v(v_points);
    const MatrixType M = mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    A.rightmultiply(M_inv);
    return XT::Common::to_string(A, precision);
  } // ... create_flux_matrix

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the moment vector is psi_vac * <b>.
  static std::string create_initial_values(const PointsVectorType& v_points = create_equidistant_points(),
                                           const RangeFieldImp psi_vac = 5e-9)
  {
    RangeType initial_vals = basisfunctions_integrated(v_points);
    initial_vals *= psi_vac;
    return XT::Common::to_string(initial_vals, precision);
  } // ... create_initial_values()

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  static void create_rhs_values(ConfigType& rhs_config, const PointsVectorType& v_points)
  {
    const RangeType integrated_basis = basisfunctions_integrated(v_points);
    const MatrixType M = mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);

    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldImp, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldImp, dimRange> b(0);
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (ii == 2) // 1 < x <= 1.5
          b[rr] = integrated_basis[rr];
        if (ii == 0 || ii == 1) // x <= 1
          A[rr][rr] = -1;
        for (size_t cc = 0; cc < dimRange; ++cc) {
          if (ii == 2 || ii == 3) // 1 <= x <= 2
            A[rr][cc] = integrated_basis[rr] * c[cc] - 3. * (rr == cc);
          else if (ii == 4 || ii == 5) // 2 <= x <= 3
            A[rr][cc] = 5 * integrated_basis[rr] * c[cc] - 10. * (rr == cc);
        } // cc
      } // rr
      rhs_config["A." + XT::Common::to_string(ii)] = XT::Common::to_string(A, precision);
      rhs_config["b." + XT::Common::to_string(ii)] = XT::Common::to_string(b, precision);
    } // ii
  } // ... create_rhs_values()

  // returns the numerator g of the left boundary value (see create_boundary_values)
  static RangeFieldImp numerator(const RangeFieldImp v)
  {
    return std::exp(-1e5 * (v - 1) * (v - 1));
  }

  // returns the denominator <g> of the left boundary value (see create_boundary_values)
  static RangeFieldImp denominator()
  {
    static constexpr double pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10) * std::erf(200 * std::sqrt(10));
  }

  // calculates integral from v_l to v_u of numerator g
  static RangeFieldImp integral_1(RangeFieldImp v_l, RangeFieldImp v_u)
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10)
           * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
  }

  static RangeFieldImp integral_2(RangeFieldImp v_l, RangeFieldImp v_u)
  {
    return integral_1(v_l, v_u) - 1. / 2e5 * (numerator(v_u) - numerator(v_l));
  }

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // Simulate with linear interpolating function.
  static std::string create_boundary_values(const PointsVectorType& v_points = create_equidistant_points(),
                                            const RangeFieldImp psi_vac = 5e-9)
  {
    const RangeType integrated_basis = basisfunctions_integrated(v_points);
    RangeType left_boundary_values(0);
    for (size_t ii = 0; ii < num_points - 1; ++ii) {
      left_boundary_values[2 * ii] = integral_1(v_points[ii], v_points[ii + 1]) / denominator();
      left_boundary_values[2 * ii + 1] = integral_2(v_points[ii], v_points[ii + 1]) / denominator();
    }

    std::string str = "[";
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0)
        str += " ";
      str += XT::Common::to_string(left_boundary_values[nn], precision) + "-("
             + XT::Common::to_string(left_boundary_values[nn] - psi_vac * integrated_basis[nn], precision)
             + ")*x[0]/3.0";
    }
    str += "]";
    return str;
  } // ... create_boundary_values()
}; // class SourceBeamPnFirstOrderDG


} // namespace Problems


template <class G, class R = double, size_t momentOrder = 5>
class SourceBeamTestCase : public Dune::GDT::Test::NonStationaryTestCase<G,
                                                                         Problems::SourceBeamPnLegendre<
                                                                             typename G::template Codim<0>::Entity,
                                                                             typename G::ctype,
                                                                             G::dimension,
                                                                             R,
                                                                             momentOrder>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;

public:
  static const size_t d = G::dimension;
  static_assert(d == 1, "Only implemented for dimension 1.");
  typedef typename Problems::SourceBeamPnLegendre<E, D, d, R, momentOrder> ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;
  using typename BaseType::LevelGridViewType;

  SourceBeamTestCase(const size_t num_refs = 1, const double divide_t_end_by = 1.0)
    : BaseType(
          divide_t_end_by, XT::Grid::make_cube_grid<GridType>(ProblemType::default_grid_config()).grid_ptr(), num_refs)
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
        << "||  Testcase: Fokker-Planck SourceBeam                                                                ||\n"
        << "|+----------------------------------------------------------------------------------------------------+|\n"
        << "||  domain = [0, 3]                                                                                   ||\n"
        << "||  time = [0, " + Dune::XT::Common::to_string(BaseType::t_end()) + "]                                ||\n"
        << "||  flux = see http://dx.doi.org/10.1137/130934210 Section 6.5                                        ||\n"
        << "||  rhs = http://dx.doi.org/10.1137/130934210 Section 6.5                                             ||\n"
        << "||  reference solution: discrete solution on finest grid                                              ||\n"
        << "|+====================================================================================================+|\n"
        << "+======================================================================================================+"
        << std::endl;
  }

private:
  const ProblemType problem_;
}; // class SourceBeamTestCase


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_SOURCEBEAM_HH
