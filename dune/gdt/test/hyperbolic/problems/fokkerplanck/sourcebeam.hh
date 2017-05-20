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


#include <dune/pdelab/common/crossproduct.hh>

#include <dune/gdt/test/instationary-eocstudy.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/math.hh>

#include "twobeams.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class TestcaseType, class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class KineticFokkerPlanck : public TwoBeamsBase<SourceBeamImp, E, D, d, R, r, rC>
{
  typedef SourceBeamBase<SourceBeamImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef TwoBeamsBase<SourceBeamImp, E, D, d, R, r, rC> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  static const size_t numPoints = 0;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return TestcaseType::static_id();
  }

  //  static ConfigType default_grid_config()
  //  {
  //    ConfigType grid_config;
  //    grid_config["type"] = "provider.cube";
  //    grid_config["lower_left"] = "[0.0]";
  //    grid_config["upper_right"] = "[3.0]";
  //    grid_config["num_elements"] = "[3000]";
  //    grid_config["overlap_size"] = "[1 1 1 1]";
  //    return grid_config;
  //  }

  using BaseType::default_boundary_info_config;

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const PointsVectorType v_points = create_equidistant_points(),
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    config.add(SourceBeamImp::create_flux_config(v_points), "flux");
    config.add(SourceBeamImp::create_rhs_config(grid_config, v_points), "rhs");
    config.add(SourceBeamImp::create_initial_value_config(grid_config, v_points, psi_vac), "initial_values");
    config.add(SourceBeamImp::create_boundary_value_config(v_points, psi_vac), "boundary_values");
    return config;
  } // ... default_config(...)

  template <class... Args>
  SourceBeamBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static bool has_non_zero_rhs()
  {
    return true;
  }

  static MatrixType mass_matrix()
  {
    return BasisFunctionType::mass_matrix();
  }

  static MatrixType mass_matrix_with_v()
  {
    return BasisFunctionType::mass_matrix_with_v();
  }

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  static std::unique_ptr<FluxType> create_flux()
  {
    MatrixType A = mass_matrix_with_v(v_points);
    const MatrixType M_inv = BasisFunctionType::mass_matrix_inverse();
    //    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    A.rightmultiply(M_inv);
    return std::make_unique<DefaultFluxType>(FluxAffineFunctionType(A));
  } // ... create_flux(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the moment vector is psi_vac * <b>.
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = SourceBeamImp::DefaultInitialValueType::static_id();
    RangeType initial_vals = SourceBeamImp::basisfunctions_integrated(v_points);
    initial_vals *= psi_vac;
    initial_value_config["values.0"] = XT::Common::to_string(initial_vals, precision);
    initial_value_config["order.0"] = "1";
    return initial_value_config;
  } // ... create_initial_values()

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  static DefaultRHSType create_rhs(const ConfigType& grid_config = default_grid_config())
  {
    const FieldVector<size_t, 3> num_elements = TestcaseType::num_elements();
    const std::vector<RangeFieldType> sigma_a = TestcaseType::sigma_a(grid_config);
    const std::vector<RangeFieldType> sigma_s = TestcaseType::sigma_s(grid_config);
    const std::vector<RangeFieldType> Q = TestcaseType::Q(grid_config);
    const size_t num_regions = std::accumulate(num_elements.begin(), num_elements.end());
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_config["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_config["upper_right"]);
    const auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = BasisFunctionType::integrated();
    const MatrixType M_inv = BasisFunctionType::mass_matrix_inverse();
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);
    MatrixType I(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType G(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        G[rr][cc] = basis_integrated[rr] * c[cc];
    const auto vol = vol();

    std::vector<RHSAffineFunctionType> affine_functions;
    for (size_t ii = 0; ii < num_regions; ++ii) {
      MatrixType G_scaled = G;
      G_scaled *= sigma_s[ii] / vol;
      MatrixType I_scaled = I;
      I_scaled *= sigma_t[ii];
      MatrixType A = G_scaled;
      A -= I_scaled;
      RangeType b = basis_integrated;
      b *= Q[ii];
      affine_functions.emplace_back(A, b);
    } // ii
    return std::make_unique<DefaultRHSType>(lower_left, upper_right, num_elements, affine_functions);
  } // ... create_rhs(...)

  // returns the numerator g of the left boundary value (see create_boundary_values)
  static RangeFieldType numerator(const RangeFieldType v)
  {
    return std::exp(-1e5 * (v - 1) * (v - 1));
  }

  // returns the denominator <g> of the left boundary value (see create_boundary_values)
  static RangeFieldType denominator()
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10) * std::erf(200 * std::sqrt(10));
  }

  // calculates integral from v_l to v_u of numerator g
  static RangeFieldType integral_1(RangeFieldType v_l, RangeFieldType v_u)
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10)
           * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
  }
}; // class SourceBeamBase


/** \see class TwoBeams in twobeams.hh */
template <class SourceBeamImp,
          class E,
          class D,
          size_t d,
          class R,
          size_t r,
          size_t rC = 1,
          class PointsVectorType = FieldVector<R, r>>
class SourceBeamBase : public TwoBeamsBase<SourceBeamImp, E, D, d, R, r, rC>
{
  typedef SourceBeamBase<SourceBeamImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef TwoBeamsBase<SourceBeamImp, E, D, d, R, r, rC> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  static const size_t numPoints = 0;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;

  static std::string static_id()
  {
    return SourceBeamImp::static_id();
  }

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

  using BaseType::default_boundary_info_config;

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const PointsVectorType v_points = create_equidistant_points(),
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    config.add(SourceBeamImp::create_flux_config(v_points), "flux");
    config.add(SourceBeamImp::create_rhs_config(grid_config, v_points), "rhs");
    config.add(SourceBeamImp::create_initial_value_config(grid_config, v_points, psi_vac), "initial_values");
    config.add(SourceBeamImp::create_boundary_value_config(v_points, psi_vac), "boundary_values");
    return config;
  } // ... default_config(...)

  template <class... Args>
  SourceBeamBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static double CFL()
  {
    return 0.4;
  }

  static double t_end()
  {
    return 4.0;
  }

  static bool has_non_zero_rhs()
  {
    return true;
  }

  static MatrixType mass_matrix(const PointsVectorType& /*v_points*/)
  {
    return MatrixType(0);
  };

  static MatrixType mass_matrix_with_v(const PointsVectorType& /*v_points*/)
  {
    return MatrixType(0);
  };

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  static ConfigType create_flux_config(const PointsVectorType& v_points = create_equidistant_points())
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A = SourceBeamImp::mass_matrix_with_v(v_points);
    const MatrixType M = SourceBeamImp::mass_matrix(v_points);
    const MatrixType M_inv = SourceBeamImp::tridiagonal_matrix_inverse(M);
    A.rightmultiply(M_inv);
    flux_config["A"] = XT::Common::to_string(A, precision);
    flux_config["b"] = Dune::XT::Common::to_string(RangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the moment vector is psi_vac * <b>.
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = SourceBeamImp::DefaultInitialValueType::static_id();
    RangeType initial_vals = SourceBeamImp::basisfunctions_integrated(v_points);
    initial_vals *= psi_vac;
    initial_value_config["values.0"] = XT::Common::to_string(initial_vals, precision);
    initial_value_config["order.0"] = "1";
    return initial_value_config;
  } // ... create_initial_values()

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 1 if x <= 2, 0 else
  // sigma_s = 0 if x <= 1, 2 if 1 < x <= 2, 10 else
  // Q = 1 if 1 <= x <= 1.5, 0 else
  static ConfigType create_rhs_config(const ConfigType& grid_config = default_grid_config(),
                                      const PointsVectorType& v_points = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[6]";
    rhs_config["name"] = static_id();
    const RangeType integrated_basis = SourceBeamImp::basisfunctions_integrated(v_points);
    const MatrixType M = SourceBeamImp::mass_matrix(v_points);
    const MatrixType M_inv = SourceBeamImp::tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);

    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldType, dimRange> b(0);
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
    return rhs_config;
  } // ... create_rhs_config(...)

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
    Dune::FieldVector<RangeFieldType, dimRange + 1> a(0), b(0), c(0), theta(0);
    Dune::FieldVector<RangeFieldType, dimRange + 2> phi(0);
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

  // returns the numerator g of the left boundary value (see create_boundary_values)
  static RangeFieldType numerator(const RangeFieldType v)
  {
    return std::exp(-1e5 * (v - 1) * (v - 1));
  }

  // returns the denominator <g> of the left boundary value (see create_boundary_values)
  static RangeFieldType denominator()
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10) * std::erf(200 * std::sqrt(10));
  }

  // calculates integral from v_l to v_u of numerator g
  static RangeFieldType integral_1(RangeFieldType v_l, RangeFieldType v_u)
  {
    static constexpr auto pi = M_PI;
    return 1 / 200. * std::sqrt(pi / 10)
           * (std::erf(100 * std::sqrt(10) * (v_u - 1)) - std::erf(100 * std::sqrt(10) * (v_l - 1)));
  }
}; // class SourceBeamBase


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_points>
class SourceBeamPnHatFunctions : public SourceBeamBase<SourceBeamPnHatFunctions<E, D, d, R, num_points>,
                                                       E,
                                                       D,
                                                       d,
                                                       R,
                                                       num_points,
                                                       1,
                                                       FieldVector<R, num_points>>
{
  typedef SourceBeamPnHatFunctions<E, D, d, R, num_points> ThisType;
  typedef SourceBeamBase<SourceBeamPnHatFunctions<E, D, d, R, num_points>,
                         E,
                         D,
                         d,
                         R,
                         num_points,
                         1,
                         FieldVector<R, num_points>>
      BaseType;

public:
  using BaseType::dimRange;
  static const size_t numPoints = num_points;
  using BaseType::precision;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::MatrixType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "SourceBeamPnHatFunctions";
  }

  template <class... Args>
  SourceBeamPnHatFunctions(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  using BaseType::create_equidistant_points;
  using BaseType::numerator;
  using BaseType::denominator;
  using BaseType::integral_1;

  static RangeType basisfunctions_integrated(const PointsVectorType& v_points = create_equidistant_points())
  {

    RangeType ret(0);
    ret[0] = v_points[1] - v_points[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      ret[ii] = v_points[ii + 1] - v_points[ii - 1];
    ret[dimRange - 1] = v_points[dimRange - 1] - v_points[dimRange - 2];
    ret *= 0.5;
    return ret;
  }

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // Simulate with linear interpolating function.
  static ConfigType create_boundary_value_config(const PointsVectorType& v_points = create_equidistant_points(),
                                                 const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    const RangeType integrated_basis = basisfunctions_integrated(v_points);
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
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config()

}; // class SourceBeamPnHatFunctions


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_points>
class SourceBeamPnFirstOrderDG : public SourceBeamBase<SourceBeamPnFirstOrderDG<E, D, d, R, num_points>,
                                                       E,
                                                       D,
                                                       d,
                                                       R,
                                                       2 * num_points - 2,
                                                       1,
                                                       FieldVector<R, num_points>>
{
  typedef SourceBeamPnFirstOrderDG<E, D, d, R, num_points> ThisType;
  typedef SourceBeamBase<SourceBeamPnFirstOrderDG<E, D, d, R, num_points>,
                         E,
                         D,
                         d,
                         R,
                         2 * num_points - 2,
                         1,
                         FieldVector<R, num_points>>
      BaseType;

public:
  using BaseType::precision;
  using BaseType::dimRange;
  static const size_t numPoints = num_points;

  using typename BaseType::ConfigType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  typedef typename Dune::FieldVector<DomainFieldType, num_points> PointsVectorType;

  static std::string static_id()
  {
    return "SourceBeamPnFirstOrderDG";
  }

  template <class... Args>
  SourceBeamPnFirstOrderDG(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  using BaseType::create_equidistant_points;

  using BaseType::numerator;
  using BaseType::denominator;
  using BaseType::integral_1;

  static RangeFieldType integral_2(RangeFieldType v_l, RangeFieldType v_u)
  {
    return integral_1(v_l, v_u) - 1. / 2e5 * (numerator(v_u) - numerator(v_l));
  }

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // Simulate with linear interpolating function.
  static ConfigType create_boundary_value_config(const PointsVectorType& v_points = create_equidistant_points(),
                                                 const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
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
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config()
}; // class SourceBeamPnFirstOrderDG


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t momentOrder>
class SourceBeamPnLegendre
    : public SourceBeamBase<SourceBeamPnLegendre<E, D, d, R, momentOrder>, E, D, d, R, momentOrder + 1, 1>
{
  typedef SourceBeamPnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef SourceBeamBase<SourceBeamPnLegendre<E, D, d, R, momentOrder>, E, D, d, R, momentOrder + 1, 1> BaseType;

public:
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::RangeFieldType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "SourceBeamPnLegendre";
  }

  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config = BaseType::default_config(grid_config, create_equidistant_points(), psi_vac);
    config.add(TwoBeamsPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder>::create_flux_config(), "flux", true);
    return config;
  } // ... default_config(...)

  template <class... Args>
  SourceBeamPnLegendre(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    RangeType ret(0);
    ret[0] = 2;
    return ret;
  }

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
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[6]";
    rhs_config["name"] = static_id();
    for (size_t ii = 0; ii < 6; ++ii) {
      MatrixType A(0);
      RangeType b(0);
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
    return rhs_config;
  } // ... create_rhs_config()

  // Boundary value of kinetic equation is \frac{exp(-10^5(v-1)^2)}{<exp(-10^5(v-1)^2)>} at x = 0 and
  // \psi_{vac} = 0.5*10^(-8) at x = 3, so n-th component of boundary value has to be
  // \frac{<base_n(v)*exp(-10^5(v-1)^2)>}{<exp(-10^5(v-1)^2)>} at x = 0 and \psi_{vac}*base_integrated_n
  // at x = 3.
  // For Legendre polynomials, this is approx. [1 0.98 0.98 ...] at x = 0 and [2*\psi_{vac} 0 0 ... ]
  // at x = 3.
  // Simulate with linear interpolating function.
  static ConfigType create_boundary_value_config(const PointsVectorType& /*v_points*/ = create_equidistant_points(),
                                                 const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    std::string str = "[" + XT::Common::to_string(left_boundary_vals[0], precision) + "-("
                      + XT::Common::to_string(left_boundary_vals[0] - 2 * psi_vac, precision) + ")*x[0]/3.0";
    for (size_t cc = 1; cc < dimRange; ++cc)
      str += " " + XT::Common::to_string(left_boundary_vals[cc], precision) + "-("
             + XT::Common::to_string(left_boundary_vals[cc], precision) + ")*x[0]/3.0";
    str += "]";
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config()

private:
  static const std::vector<R> left_boundary_vals;
}; // class SourceBeamPnLegendre<...>

template <class E, class D, size_t d, class R, size_t momentOrder>
const std::vector<R> SourceBeamPnLegendre<E, D, d, R, momentOrder>::left_boundary_vals = {
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
template <class E, class D, size_t d, class R, size_t momentOrder>
class SourceBeamPnLegendreLaplaceBeltrami
    : public SourceBeamBase<SourceBeamPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder>,
                            E,
                            D,
                            d,
                            R,
                            momentOrder + 1,
                            1>
{
  typedef SourceBeamPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder> ThisType;
  typedef SourceBeamBase<SourceBeamPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder>, E, D, d, R, momentOrder + 1, 1>
      BaseType;

public:
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "SourceBeamPnLegendreLaplaceBeltrami";
  }

public:
  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const RangeFieldType psi_vac = 1e-4)
  {
    ConfigType config = BaseType::default_config(grid_config, create_equidistant_points(), psi_vac);
    config.add(TwoBeamsPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder>::create_flux_config(), "flux", true);
    return config;
  } // ... default_config(...)

  template <class... Args>
  SourceBeamPnLegendreLaplaceBeltrami(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnLegendre<E, D, d, R, momentOrder>::basisfunctions_integrated(v_points);
  }

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
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[6]";
    rhs_config["name"] = DefaultRHSType::static_id();
    for (size_t ii = 0; ii < 6; ++ii) {
      Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> A(0);
      Dune::FieldVector<RangeFieldType, dimRange> b(0);
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
    return rhs_config;
  } // ... create_rhs_config()

  // Boundary value of kinetic equation is delta(v-1) at x = 0 and psi_vac at x = 3,
  // so n-th component of boundary value has to be 0.5*\phi_n(1) at x = 0 and 0.5*psi_vac*delta(n) at x = 3.
  // For Legendre polynomials, this is [0.5 0.5 0.5 ...] at x = 0 and [2*psi_vac 0 0 ... ] at x = 3.
  // Model with linear interpolating function.
  static ConfigType create_boundary_value_config(const PointsVectorType& /*v_points*/ = create_equidistant_points(),
                                                 const RangeFieldType psi_vac = 1e-4)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    std::string str = "[0.5-" + XT::Common::to_string(0.5 - 2 * psi_vac, precision) + "*x[0]/3.0";
    for (size_t cc = 1; cc < dimRange; ++cc)
      str += " 0.5-0.5*x[0]/3.0";
    str += "]";
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_values()
}; // class SourceBeamPnLegendreLaplaceBeltrami


} // namespace Problems


template <class G, class R = double, size_t momentOrder = 5>
class SourceBeamTestCase : public Dune::GDT::Test::NonStationaryTestCase<G,
                                                                         Problems::SourceBeamPnLegendreLaplaceBeltrami<
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
  typedef typename Problems::SourceBeamPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder> ProblemType;
  static const size_t dimRange = ProblemType::dimRange;
  static const size_t dimRangeCols = 1;

private:
  typedef typename Dune::GDT::Test::NonStationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::SolutionType;

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
