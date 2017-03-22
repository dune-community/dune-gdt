// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_CHECKERBOARD3D_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_CHECKERBOARD3D_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include <dune/xt/functions/global.hh>

#include "pointsource.hh"


namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

/** \see class TwoBeams in twobeams.hh */
template <class CheckerboardImp,
          class E,
          class D,
          size_t d,
          class R,
          size_t r,
          size_t rC = 1,
          class PointsVectorType = FieldVector<R, r>>
class CheckerboardBase : public PlaneSourceBase<CheckerboardImp, E, D, d, R, r, rC, PointsVectorType>
{
  typedef CheckerboardBase<CheckerboardImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef PlaneSourceBase<CheckerboardImp, E, D, d, R, r, rC, PointsVectorType> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[0 0 0]";
    grid_config["upper_right"] = "[7 7 7]";
    grid_config["num_elements"] = "[70 70 70]";
    grid_config["overlap_size"] = "[1 1 1 1]";
    return grid_config;
  }

  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;
  using BaseType::tridiagonal_matrix_inverse;

  template <class... Args>
  CheckerboardBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static double CFL()
  {
    return 0.4;
  }

  static double t_end()
  {
    return 3.2;
  }

  // Initial value of the kinetic equation is psi_vac.
  // Thus the initial value for the n-th moment is base_integrated_n * psi_vac
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = CheckerboardImp::DefaultInitialValueType::static_id();
    const auto basis_integrated = CheckerboardImp::basisfunctions_integrated(v_points);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr] * psi_vac, precision);
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "1";
    return initial_value_config;
  } // ... create_initial_value_config(...)

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& v_points = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = CheckerboardImp::DefaultRHSType::static_id();
    const RangeType integrated_basis = CheckerboardImp::basisfunctions_integrated(v_points);
    const MatrixType M = CheckerboardImp::mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);
    Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> A(0);
    Dune::FieldVector<RangeFieldType, dimRange> b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        A[rr][cc] = 0.5 * integrated_basis[rr] * c[cc] - (rr == cc);
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()
}; // class PointSourceBase


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t momentOrder>
class CheckerboardPnLegendre : public PointSourceBase<CheckerboardPnLegendre<E, D, d, R, momentOrder>,
                                                      E,
                                                      D,
                                                      d,
                                                      R,
                                                      (momentOrder + 1) * (momentOrder + 1),
                                                      1>
{
  typedef CheckerboardPnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef PointSourceBase<CheckerboardPnLegendre<E, D, d, R, momentOrder>,
                          E,
                          D,
                          d,
                          R,
                          (momentOrder + 1) * (momentOrder + 1),
                          1>
      BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::FluxRangeType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "CheckerboardPnLegendre";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;

  template <class... Args>
  CheckerboardPnLegendre(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // flux matrix A_i,nm = <Omega_i h_n h_m>
  static ConfigType create_flux_config(const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A_0(0), A_1(0), A_2(0);
    const auto quadrature = get_lebedev_quadrature(131);
    for (int nn = 0; nn < dimRange; ++nn) {
      for (int mm = 0; mm < dimRange; ++mm) {
        auto l_and_m_nn = get_l_and_m(nn);
        auto l_and_m_mm = get_l_and_m(mm);
        for (const auto& quad_point : quadrature) {
          const auto point = quad_point.position();
          const auto mu = point[0];
          const auto phi = point[1];
          const auto weight = quad_point.weight();
          A_0[nn][mm] += std::sqrt(1. - mu * mu) * std::cos(phi)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
          A_1[nn][mm] += std::sqrt(1. - mu * mu) * std::sin(phi)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
          A_2[nn][mm] += mu * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
        }
      }
    }
    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    RangeType ret(0);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  // n-th component of RHS is -sigma_t u_n + (sigma_s u_0 + sqrt(4*pi)*Q) delta(n)
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // Thus, rhs[n] = (delta(n)-1)u[n]
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = DefaultRHSType::static_id();
    MatrixType A(0);
    RangeType b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      A[rr][rr] = (rr == 0) - 1;
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()
}; // class PointSourcePnLegendre


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_vertices>
class CheckerboardPnHatFunctions
    : public CheckerboardBase<CheckerboardPnHatFunctions<E, D, d, R, num_vertices>, E, D, d, R, num_vertices, 1>
{
  typedef CheckerboardPnHatFunctions<E, D, d, R, num_vertices> ThisType;
  typedef CheckerboardBase<CheckerboardPnHatFunctions<E, D, d, R, num_vertices>, E, D, d, R, num_vertices, 1> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RHSType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultRHSType;
  typedef typename XT::Functions::GlobalLambdaFunction<E, D, d, R, dimRange> InitialValueFunctionType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<InitialValueFunctionType, E, D, d, R, dimRange>
      DefaultInitialValueType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::DomainType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "CheckerboardPnHatFunctions";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;

  static std::unique_ptr<ThisType> create(const Dune::QuadratureRule<double, 3>& quadrature,
                                          const SphericalTriangulation<double>& poly,
                                          const ConfigType config = default_config())
  {
    const std::shared_ptr<const FluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const RHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const InitialValueType> initial_values(DefaultInitialValueType::create(
        config.sub("initial_values"), "", create_initial_value_lambda(quadrature, poly)));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const BoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const ConfigType grid_config,
                                   const Dune::QuadratureRule<double, 3>& quadrature,
                                   const SphericalTriangulation<double>& poly,
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    config.add(create_flux_config(quadrature, poly), "flux");
    config.add(create_rhs_config(grid_config, quadrature, poly), "rhs");
    config.add(create_initial_value_config(grid_config, quadrature, poly, psi_vac), "initial_values");
    config.add(create_boundary_value_config(quadrature, poly, psi_vac), "boundary_values");
    return config;
  } // ... default_config(...)

  template <class... Args>
  CheckerboardPnHatFunctions(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // flux matrix A_i,nm = <Omega_i h_n h_m>
  static ConfigType create_flux_config(const Dune::QuadratureRule<double, 3>& quadrature,
                                       const SphericalTriangulation<double>& poly)
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A_0(0), A_1(0), A_2(0);
    for (const auto& quad_point : quadrature) {
      const auto& point = quad_point.position();
      const auto basis_evaluated = evaluate_spherical_barycentric_coordinates<RangeType, DomainType>(point, poly);
      const auto& weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn) {
        for (size_t mm = 0; mm < dimRange; ++mm) {
          A_0[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[0] * weight;
          A_1[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[1] * weight;
          A_2[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[2] * weight;
        } // mm
      } // nn
    } // quadrature
    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const Dune::QuadratureRule<double, 3>& quadrature,
                                             const SphericalTriangulation<double>& poly)
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature) {
      const auto point = quad_point.position();
      const auto basis_evaluated = evaluate_spherical_barycentric_coordinates<RangeType, DomainType>(point, poly);
      const auto weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        ret[nn] += basis_evaluated[nn] * weight;
    } // quadrature
    return ret;
  }

  static bool is_absorbing(size_t plane, size_t row, size_t col)
  {
    if (plane == 0 || plane == 6)
      return false;
    if (plane == 2 || plane == 4)
      return !((row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
               || ((row == 3 || row == 5) && (col == 1 || col == 5)))
             && (row != 0 && row != 6) && (col != 0 && col != 6) && !(row == 5 && col == 3) && !(row == 3 && col == 3);
    if (plane == 3 || plane == 1 || plane == 5)
      return (row == 1 && col % 2 == 1) || ((row == 2 || row == 4) && (col == 2 || col == 4))
             || ((row == 3 || row == 5) && (col == 1 || col == 5))
             || (plane != 3 && col == 3 && (row == 3 && row == 5));
  }

  // n-th component of RHS is -sigma_t u_n + sigma_s/(4 pi) <psi><b_n> + Q<b_n>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0 in scattering regions
  // sigma_a = 10, sigma_s = 0, Q = 0 in absorbing regions
  // sigma_a = 0, sigma_s = 1, Q = 1 for center cube
  // As <psi> = sum_i u_i, the rhs becomes
  // - sigma_t u_n + sigma_s/(4 pi) sum_i u_i <b_n> + Q<b_n>
  static ConfigType create_rhs_config(const ConfigType grid_config,
                                      const Dune::QuadratureRule<double, 3>& quadrature,
                                      const SphericalTriangulation<double>& poly)
  {
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[7 7 7]";
    rhs_config["name"] = DefaultRHSType::static_id();

    RangeType Q;
    RangeFieldType sigma_s, sigma_t;

    for (size_t plane = 0; plane < 7; ++plane) {
      for (size_t row = 0; row < 7; ++row) {
        for (size_t col = 0; col < 7; ++col) {
          if (plane == 3 && row == 3 && col == 3) { // center
            Q = basis_integrated;
            sigma_s = 1;
            sigma_t = 1;
          } else if (is_absorbing(plane, row, col)) { // absorbing regions
            Q *= 0;
            sigma_s = 0;
            sigma_t = 10;
          } else { // scattering regions (without center)
            Q *= 0;
            sigma_s = 1;
            sigma_t = 1;
          }

          MatrixType A(1);
          A *= sigma_s / (4 * M_PI);
          for (size_t nn = 0; nn < dimRange; ++nn) {
            A[nn] *= basis_integrated[nn];
            A[nn][nn] -= sigma_t;
          }
          size_t number = 49 * plane + 7 * row + col;
          rhs_config["A." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(A, precision);
          rhs_config["b." + Dune::XT::Common::to_string(number)] = Dune::XT::Common::to_string(Q, precision);
        } // col
      } // row
    } // plane

    return rhs_config;
  } // ... create_rhs_config()

  // Initial value of the kinetic equation is psi_vac
  static ConfigType create_initial_value_config(const ConfigType& grid_config,
                                                const Dune::QuadratureRule<double, 3>& quadrature,
                                                const SphericalTriangulation<double>& poly,
                                                const double psi_vac = 5e-9)
  {
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = DefaultInitialValueType::static_id();
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr] * psi_vac, precision);
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "1";
    return initial_value_config;
  } // ... create_initial_value_config(...)

  // Initial value of the kinetic equation is psi_vac
  // Thus the initial value for the n-th moment is base_integrated_n * psi_vac
  static std::vector<std::function<RangeType(DomainType)>>
  create_initial_value_lambda(const Dune::QuadratureRule<double, 3>& quadrature,
                              const SphericalTriangulation<double>& poly,
                              const double psi_vac = 5e-9)
  {
    std::vector<std::function<RangeType(DomainType)>> ret;
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ret.push_back([basis_integrated, psi_vac](const DomainType& x) {
      RangeType result = basis_integrated;
      result *= psi_vac;
      return result;
    });
    return ret;
  } // ... create_initial_value_lambda(...)

  // Boundary value of kinetic equation is psi_vac at all boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at all boundaries.
  // Modell with constant function.
  static ConfigType create_boundary_value_config(const Dune::QuadratureRule<double, 3>& quadrature,
                                                 const SphericalTriangulation<double>& poly,
                                                 const double psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    const RangeType integrated_basis = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(integrated_basis[rr] * psi_vac, precision);
    }
    str += "]";
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config(...)

}; // class CheckerboardPnHatFunctions


///** \see class TwoBeams in twobeams.hh */
// template <class E, class D, size_t d, class R, size_t num_faces>
// class PointSourcePnPartialMoments
//    : public PointSourceBase<PointSourcePnPartialMoments<E, D, d, R, num_faces>, E, D, d, R, 4 * num_faces, 1>
//{
//  typedef PointSourcePnPartialMoments<E, D, d, R, num_faces> ThisType;
//  typedef PointSourceBase<PointSourcePnPartialMoments<E, D, d, R, num_faces>, E, D, d, R, 4 * num_faces, 1> BaseType;

// public:
//  using BaseType::dimDomain;
//  using BaseType::dimRange;
//  using BaseType::precision;
//  using typename BaseType::FluxRangeType;
//  using typename BaseType::RHSType;
//  using typename BaseType::FluxType;
//  using typename BaseType::InitialValueType;
//  using typename BaseType::BoundaryValueType;
//  using typename BaseType::DefaultFluxType;
//  using typename BaseType::DefaultRHSType;
//  typedef typename XT::Functions::GlobalLambdaFunction<E, D, d, R, dimRange> InitialValueFunctionType;
//  typedef typename XT::Functions::FunctionCheckerboardFunction<InitialValueFunctionType, E, D, d, R, dimRange>
//      DefaultInitialValueType;
//  using typename BaseType::DefaultBoundaryValueType;
//  using typename BaseType::ConfigType;
//  using typename BaseType::MatrixType;
//  using typename BaseType::RangeFieldType;
//  using typename BaseType::RangeType;
//  using typename BaseType::DomainType;
//  typedef RangeType PointsVectorType;

//  static std::string static_id()
//  {
//    return "PointSourcePnPartialBasis";
//  }

//  using BaseType::create_equidistant_points;
//  using BaseType::default_grid_config;
//  using BaseType::default_boundary_info_config;

//  static std::unique_ptr<ThisType> create(const Dune::QuadratureRule<double, 3>& quadrature,
//                                          const SphericalTriangulation<double>& poly,
//                                          const ConfigType config = default_config())
//  {
//    const std::shared_ptr<const FluxType> flux(DefaultFluxType::create(config.sub("flux")));
//    const std::shared_ptr<const RHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
//    const std::shared_ptr<const InitialValueType> initial_values(DefaultInitialValueType::create(
//        config.sub("initial_values"), "", create_initial_value_lambda(quadrature, poly)));
//    const ConfigType grid_config = config.sub("grid");
//    const ConfigType boundary_info = config.sub("boundary_info");
//    const std::shared_ptr<const BoundaryValueType> boundary_values(
//        DefaultBoundaryValueType::create(config.sub("boundary_values")));
//    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
//  } // ... create(...)

//  static ConfigType default_config(const ConfigType grid_config,
//                                   const Dune::QuadratureRule<double, 3>& quadrature,
//                                   const SphericalTriangulation<double>& poly,
//                                   const RangeFieldType psi_vac = 5e-9)
//  {
//    ConfigType config;
//    config.add(grid_config, "grid");
//    config.add(default_boundary_info_config(), "boundary_info");
//    config.add(create_flux_config(quadrature, poly), "flux");
//    config.add(create_rhs_config(grid_config, quadrature, poly), "rhs");
//    config.add(create_initial_value_config(grid_config, quadrature, poly, psi_vac), "initial_values");
//    config.add(create_boundary_value_config(quadrature, poly, psi_vac), "boundary_values");
//    return config;
//  } // ... default_config(...)

//  template <class... Args>
//  PointSourcePnPartialMoments(Args&&... args)
//    : BaseType(std::forward<Args>(args)...)
//  {
//  }

//  // flux matrix A_i,nm = <Omega_i h_n h_m>
//  static ConfigType create_flux_config(const Dune::QuadratureRule<double, 3>& quadrature,
//                                       const SphericalTriangulation<double>& poly)
//  {
//    ConfigType flux_config;
//    flux_config["type"] = DefaultFluxType::static_id();
//    MatrixType A_0(0), A_1(0), A_2(0);
//    for (const auto& quad_point : quadrature) {
//      const auto point = quad_point.position();
//      const auto basis_evaluated = evaluate_linear_partial_basis<RangeType, DomainType>(point, poly);
//      const auto weight = quad_point.weight();
//      for (size_t nn = 0; nn < dimRange; ++nn) {
//        for (size_t mm = 0; mm < dimRange; ++mm) {
//          A_0[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[0] * weight;
//          A_1[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[1] * weight;
//          A_2[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[2] * weight;
//        } // mm
//      } // nn
//    } // quadrature
//    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
//    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
//    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
//    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
//    return flux_config;
//  } // ... create_flux_config(...)

//  // returns <b>, where b is the basis functions vector
//  static RangeType basisfunctions_integrated(const Dune::QuadratureRule<double, 3>& quadrature,
//                                             const SphericalTriangulation<double>& poly)
//  {
//    RangeType ret(0);
//    for (const auto& quad_point : quadrature) {
//      const auto& point = quad_point.position();
//      const auto basis_evaluated = evaluate_linear_partial_basis<RangeType, DomainType>(point, poly);
//      const auto& weight = quad_point.weight();
//      for (size_t nn = 0; nn < dimRange; ++nn)
//        ret[nn] += basis_evaluated[nn] * weight;
//    } // quadrature
//    return ret;
//  }

//  // n-th component of RHS is -sigma_t u_n + sigma_s/(4 pi) <psi><b_n> + Q<b_n>
//  // For this test case (sigma_t = sigma_s + sigma_a),
//  // sigma_a = 0, sigma_s = 1, Q = 0
//  // As <psi> = sum_i u_i, the rhs becomes
//  // -u_n + 1/(4 pi) sum_i u_i <b_n>
//  static ConfigType create_rhs_config(const ConfigType grid_config,
//                                      const Dune::QuadratureRule<double, 3>& quadrature,
//                                      const SphericalTriangulation<double>& poly)
//  {
//    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
//    ConfigType rhs_config;
//    rhs_config["lower_left"] = grid_config["lower_left"];
//    rhs_config["upper_right"] = grid_config["upper_right"];
//    rhs_config["num_elements"] = "[1 1 1]";
//    rhs_config["name"] = DefaultRHSType::static_id();
//    MatrixType A(1);
//    A *= 1. / (4 * M_PI);
//    RangeType b(0);
//    for (size_t nn = 0; nn < dimRange; ++nn) {
//      A[nn] *= basis_integrated[nn];
//      A[nn][nn] -= 1.;
//    }
//    rhs_config["A.0"] = XT::Common::to_string(A);
//    rhs_config["b.0"] = XT::Common::to_string(b);
//    return rhs_config;
//  } // ... create_rhs_config()

//  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
//  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
//  // exp(-|x|^2/(2*sigma^2))).
//  static ConfigType create_initial_value_config(const ConfigType& grid_config,
//                                                const Dune::QuadratureRule<double, 3>& quadrature,
//                                                const SphericalTriangulation<double>& poly,
//                                                const double psi_vac = 5e-9)
//  {
//    static const double sigma = 0.03;
//    ConfigType initial_value_config;
//    initial_value_config["lower_left"] = grid_config["lower_left"];
//    initial_value_config["upper_right"] = grid_config["upper_right"];
//    initial_value_config["num_elements"] = "[1 1 1]";
//    initial_value_config["variable"] = "x";
//    initial_value_config["name"] = DefaultInitialValueType::static_id();
//    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
//    std::string str = "[";
//    for (size_t rr = 0; rr < dimRange; ++rr) {
//      if (rr > 0)
//        str += " ";
//      str += XT::Common::to_string(basis_integrated[rr], precision) + "*(" + XT::Common::to_string(psi_vac, precision)
//             + "+" + XT::Common::to_string(1. / (8. * M_PI * sigma * sigma), precision)
//             + "*exp((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/(" + XT::Common::to_string(-2. * sigma * sigma, precision) +
//             ")))";
//    } // rr
//    str += "]";
//    initial_value_config["values.0"] = str;
//    initial_value_config["order.0"] = "20";
//    return initial_value_config;
//  } // ... create_initial_value_config(...)

//  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
//  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
//  // exp(-|x|^2/(2*sigma^2))).
//  static std::vector<std::function<RangeType(DomainType)>>
//  create_initial_value_lambda(const Dune::QuadratureRule<double, 3>& quadrature,
//                              const SphericalTriangulation<double>& poly,
//                              const double psi_vac = 5e-9)
//  {
//    std::vector<std::function<RangeType(DomainType)>> ret;
//    static const double sigma = 0.03;
//    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
//    ret.push_back([basis_integrated, psi_vac](const DomainType& x) {
//      RangeType result = basis_integrated;
//      result *= psi_vac + 1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm() / (2. * sigma * sigma));
//      return result;
//    });
//    return ret;
//  } // ... create_initial_value_lambda(...)

//  // Boundary value of kinetic equation is psi_vac at both boundaries
//  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
//  // Modell with constant function.
//  static ConfigType create_boundary_value_config(const Dune::QuadratureRule<double, 3>& quadrature,
//                                                 const SphericalTriangulation<double>& poly,
//                                                 const double psi_vac = 5e-9)
//  {
//    ConfigType boundary_value_config;
//    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
//    boundary_value_config["variable"] = "x";
//    boundary_value_config["order"] = "1";
//    const RangeType integrated_basis = basisfunctions_integrated(quadrature, poly);
//    std::string str = "[";
//    for (size_t rr = 0; rr < dimRange; ++rr) {
//      if (rr > 0)
//        str += " ";
//      str += XT::Common::to_string(integrated_basis[rr] * psi_vac, precision);
//    }
//    str += "]";
//    boundary_value_config["expression"] = str;
//    return boundary_value_config;
//  } // ... create_boundary_value_config(...)

//}; // class PointSourcePnPartialMoments


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_CHECKERBOARD3D_HH
