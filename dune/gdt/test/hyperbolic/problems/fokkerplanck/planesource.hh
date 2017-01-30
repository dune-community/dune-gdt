// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include "sourcebeam.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


/** \see class TwoBeams in twobeams.hh */
template <class PlaneSourceImp,
          class E,
          class D,
          size_t d,
          class R,
          size_t r,
          size_t rC = 1,
          class PointsVectorType = FieldVector<R, r>>
class PlaneSourceBase : public SourceBeamBase<PlaneSourceImp, E, D, d, R, r, rC, PointsVectorType>
{
  typedef PlaneSourceBase<PlaneSourceImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef SourceBeamBase<PlaneSourceImp, E, D, d, R, r, rC, PointsVectorType> BaseType;

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
    grid_config["lower_left"] = "[-1.2]";
    grid_config["upper_right"] = "[1.2]";
    grid_config["num_elements"] = "[2400]";
    grid_config["overlap_size"] = "[1 1 1 1]";
    return grid_config;
  }

  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;
  using BaseType::tridiagonal_matrix_inverse;

  template <class... Args>
  PlaneSourceBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static double CFL()
  {
    return 0.4;
  }

  static double t_end()
  {
    return 1.0;
  }

  // Initial value of the kinetic equation is psi_vac + delta(x).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + delta(x))
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = grid_config["num_elements"];
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = PlaneSourceImp::DefaultInitialValueType::static_id();
    const auto basis_integrated = PlaneSourceImp::basisfunctions_integrated(v_points);
    const size_t num_elements = XT::Common::from_string<std::vector<size_t>>(grid_config["num_elements"])[0];
    const RangeFieldType lower_left = XT::Common::from_string<std::vector<double>>(grid_config["lower_left"])[0];
    const RangeFieldType upper_right = XT::Common::from_string<std::vector<double>>(grid_config["upper_right"])[0];
    const RangeFieldType vol_entity = (upper_right - lower_left) / num_elements;
    if (num_elements % 2)
      DUNE_THROW(Dune::NotImplemented, "An even number of grid cells is needed for this test!");
    for (size_t ii = 0; ii < num_elements; ++ii) {
      std::string str = "[";
      for (size_t rr = 0; rr < dimRange; ++rr) {
        if (rr > 0)
          str += " ";
        // approximate delta function by constant value of 1/(2*vol_entity) on cells on both side of 0.
        if (ii == num_elements / 2 || ii == num_elements / 2 - 1)
          str += XT::Common::to_string(basis_integrated[rr] * (psi_vac + 1. / (2. * vol_entity)), precision);
        else
          str += XT::Common::to_string(basis_integrated[rr] * psi_vac, precision);
      } // rr
      str += "]";
      initial_value_config["values." + XT::Common::to_string(ii)] = str;
      initial_value_config["order." + XT::Common::to_string(ii)] = 1;
    } // ii
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
    rhs_config["num_elements"] = "[1]";
    rhs_config["name"] = PlaneSourceImp::DefaultRHSType::static_id();
    const RangeType integrated_basis = PlaneSourceImp::basisfunctions_integrated(v_points);
    const MatrixType M = PlaneSourceImp::mass_matrix(v_points);
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

  // Boundary value of kinetic equation is psi_vac at both boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
  // Modell with constant function.
  static ConfigType create_boundary_value_config(const RangeType& v_points = create_equidistant_points(),
                                                 const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = PlaneSourceImp::DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    const RangeType integrated_basis = PlaneSourceImp::basisfunctions_integrated(v_points);
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
}; // class PlaneSourceBase


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t momentOrder>
class PlaneSourcePnLegendre
    : public PlaneSourceBase<PlaneSourcePnLegendre<E, D, d, R, momentOrder>, E, D, d, R, momentOrder + 1, 1>
{
  typedef PlaneSourcePnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef PlaneSourceBase<PlaneSourcePnLegendre<E, D, d, R, momentOrder>, E, D, d, R, momentOrder + 1, 1> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "PlaneSourcePnLegendre";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;

  static ConfigType default_config(const ConfigType grid_config = default_grid_config(),
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config = BaseType::default_config(grid_config, create_equidistant_points(), psi_vac);
    config.add(TwoBeamsPnLegendreLaplaceBeltrami<E, D, d, R, momentOrder>::create_flux_config(), "flux", true);
    return config;
  } // ... default_config(...)

  template <class... Args>
  PlaneSourcePnLegendre(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnLegendre<E, D, d, R, momentOrder>::basisfunctions_integrated(v_points);
  }

  // n-th component of RHS is -sigma_t u_n + (sigma_s u_0 + 2Q) delta(n)
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // Thus, rhs[n] = (delta(n)-1)u[n]
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1]";
    rhs_config["name"] = DefaultRHSType::static_id();
    MatrixType A(0);
    RangeType b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      A[rr][rr] = (rr == 0) - 1;
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()
}; // class PlaneSourcePnLegendre

/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_points>
class PlaneSourcePnHatFunctions : public PlaneSourceBase<PlaneSourcePnHatFunctions<E, D, d, R, num_points>,
                                                         E,
                                                         D,
                                                         d,
                                                         R,
                                                         num_points,
                                                         1,
                                                         FieldVector<R, num_points>>
{
  typedef PlaneSourcePnHatFunctions<E, D, d, R, num_points> ThisType;
  typedef PlaneSourceBase<PlaneSourcePnHatFunctions<E, D, d, R, num_points>,
                          E,
                          D,
                          d,
                          R,
                          num_points,
                          1,
                          FieldVector<R, num_points>>
      BaseType;

public:
  static const size_t numPoints = num_points;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  typedef FieldVector<RangeFieldType, num_points> PointsVectorType;

  static std::string static_id()
  {
    return "PlaneSourcePnHatFunctions";
  }

  template <class... Args>
  PlaneSourcePnHatFunctions(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  using BaseType::create_equidistant_points;

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnHatFunctions<E, D, d, R, num_points>::basisfunctions_integrated(v_points);
  }

  static MatrixType mass_matrix(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnHatFunctions<E, D, d, R, num_points>::mass_matrix(v_points);
  }

  static MatrixType mass_matrix_with_v(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnHatFunctions<E, D, d, R, num_points>::mass_matrix_with_v(v_points);
  }
}; // class PlaneSourcePnHatFunctions


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_points>
class PlaneSourcePnFirstOrderDG : public PlaneSourceBase<PlaneSourcePnFirstOrderDG<E, D, d, R, num_points>,
                                                         E,
                                                         D,
                                                         d,
                                                         R,
                                                         2 * num_points - 2,
                                                         1,
                                                         FieldVector<R, num_points>>
{
  typedef PlaneSourcePnFirstOrderDG<E, D, d, R, num_points> ThisType;
  typedef PlaneSourceBase<PlaneSourcePnFirstOrderDG<E, D, d, R, num_points>,
                          E,
                          D,
                          d,
                          R,
                          2 * num_points - 2,
                          1,
                          FieldVector<R, num_points>>
      BaseType;

public:
  static const size_t numPoints = num_points;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  typedef FieldVector<RangeFieldType, num_points> PointsVectorType;

  static std::string static_id()
  {
    return "PlaneSourcePnFirstOrderDG";
  }

  template <class... Args>
  PlaneSourcePnFirstOrderDG(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  using BaseType::create_equidistant_points;

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnFirstOrderDG<E, D, d, R, num_points>::basisfunctions_integrated(v_points);
  }

  static MatrixType mass_matrix(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnFirstOrderDG<E, D, d, R, num_points>::mass_matrix(v_points);
  }

  static MatrixType mass_matrix_with_v(const PointsVectorType& v_points = create_equidistant_points())
  {
    return SourceBeamPnFirstOrderDG<E, D, d, R, num_points>::mass_matrix_with_v(v_points);
  }
}; // class PlaneSourcePnHatFunctions


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_PLANESOURCE_HH
