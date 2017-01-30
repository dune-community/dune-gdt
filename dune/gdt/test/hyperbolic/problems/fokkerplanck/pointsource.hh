// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "planesource.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

inline std::string trim(const std::string& s)
{
  auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
  auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
  return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

static const std::vector<int> allowed_degrees = {6,    14,   26,   38,   50,   74,   86,   110,  146,  170,  194,
                                                 230,  266,  302,  350,  434,  590,  770,  974,  1202, 1454, 1730,
                                                 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
static const std::vector<int> allowed_orders = {3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
                                                41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};

Dune::QuadratureRule<double, 2> get_lebedev_quadrature(size_t requested_order)
{
  size_t index = -1;
  for (int ii = 0; ii < allowed_orders.size(); ++ii) {
    if (allowed_orders[ii] >= requested_order) {
      index = ii;
      break;
    }
  }
  int order = allowed_orders[index];
  //  int degree = allowed_degrees[index];
  char orderstring[4];
  sprintf(orderstring, "%03d", order);
  std::string filename =
      std::string("/home/tobias/Software/dune-gdt-super-2.5/dune-gdt/dune/gdt/LebedevTables/lebedev_") + orderstring
      + ".txt";

  Dune::QuadratureRule<double, 2> quadrature_rule;
  std::string current_line;
  std::ifstream quadrature_file(filename);
  while (getline(quadrature_file, current_line)) {
    current_line = trim(current_line);
    auto quadrature_values =
        XT::Common::tokenize(current_line, " ", boost::algorithm::token_compress_mode_type::token_compress_on);
    std::cout << XT::Common::to_string(quadrature_values) << std::endl;
    const double phi = XT::Common::from_string<double>(quadrature_values[0]) / 360. * 2 * M_PI;
    const double theta = XT::Common::from_string<double>(quadrature_values[1]) / 360. * 2 * M_PI;
    const double mu = std::cos(theta);
    const double weight = 4 * M_PI * XT::Common::from_string<double>(quadrature_values[2]);
    Dune::QuadraturePoint<double, 2> quadrature_point(FieldVector<double, 2>{mu, phi}, weight);
    quadrature_rule.push_back(quadrature_point);
  }
  return quadrature_rule;
}


/** \see class TwoBeams in twobeams.hh */
template <class PointSourceImp,
          class E,
          class D,
          size_t d,
          class R,
          size_t r,
          size_t rC = 1,
          class PointsVectorType = FieldVector<R, r>>
class PointSourceBase : public PlaneSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType>
{
  typedef PointSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef PlaneSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType> BaseType;

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
    grid_config["lower_left"] = "[-0.5 -0.5 -0.5]";
    grid_config["upper_right"] = "[0.5 0.5 0.5]";
    grid_config["num_elements"] = "[100 100 100]";
    grid_config["overlap_size"] = "[1 1 1 1]";
    return grid_config;
  }

  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;
  using BaseType::tridiagonal_matrix_inverse;

  template <class... Args>
  PointSourceBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static double CFL()
  {
    return 0.4;
  }

  static double t_end()
  {
    return 0.45;
  }

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    static const double sigma = 0.03;
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = PointSourceImp::DefaultInitialValueType::static_id();
    const auto basis_integrated = PointSourceImp::basisfunctions_integrated(v_points);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr], precision) + "*(" + XT::Common::to_string(psi_vac, precision)
             + "+" + XT::Common::to_string(1. / (8. * M_PI * sigma * sigma), precision)
             + "*exp((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/(" + XT::Common::to_string(-2. * sigma * sigma, precision) + ")))";
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "20";
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
    rhs_config["name"] = PointSourceImp::DefaultRHSType::static_id();
    const RangeType integrated_basis = PointSourceImp::basisfunctions_integrated(v_points);
    const MatrixType M = PointSourceImp::mass_matrix(v_points);
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
}; // class PlaneSourceBase


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t momentOrder>
class PointSourcePnLegendre : public PointSourceBase<PointSourcePnLegendre<E, D, d, R, momentOrder>,
                                                     E,
                                                     D,
                                                     d,
                                                     R,
                                                     (momentOrder + 1) * (momentOrder + 1),
                                                     1>
{
  typedef PointSourcePnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef PointSourceBase<PointSourcePnLegendre<E, D, d, R, momentOrder>,
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
    return "PointSourcePnLegendre";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;

  template <class... Args>
  PointSourcePnLegendre(Args&&... args)
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


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
