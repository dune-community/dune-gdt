// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_TOOLS_EULER_HH
#define DUNE_GDT_TOOLS_EULER_HH

#include <cmath>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/base/sliced.hh>
#include <dune/xt/functions/base/transformed.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


/**
 * \brief A collection of useful tools in the context of the Euler equations.
 *
 * In paritcular, provides
 *
 * conversion between primitive (physical) variables
 *   - density rho: R^d -> R
 *   - velocity  v: R^d -> R^d
 *   - pressure  p: R^d -> R
 * and conservative variables w = (rho, rho*v[1], ..., rho*v[d], E): R^d -> R^m, with m = d + 2 and
 *   - total energy E: R^d -> R;
 *
 * evaluation of the euler fluxes f_1, ..., f_d: R^m -> R^m and their respective jacobians A_1, ..., A_d;
 *
 * the eigendecomposition T^{-1} P T = diagonal(\lambda_1, ..., \lambda_m) with
 *   - P = A * n (n being a unit normal)
 *   - the eigenvalues \lambda_i of P
 *   - the (right) eigenvectors T of P
 *   - the (left eigenvectors) inverse of the (right) eigenvectors T^{-1} of P;
 *
 * some quantities of interest
 *   - speed of sound
 *   - mach number
 *   - enthalpy;
 *
 * a means to visualize a function representing the conservative variables w.
 *
 * \sa [Dolejsi, Feistauer, 2016, pp. 404 and following]
 * \sa [Kröner, 1997, p. 387] for 1d and 2d
 */
template <size_t d, class R = double>
class EulerTools
{
  static_assert(1 <= d && d <= 3, "");

public:
  static const constexpr size_t m = d + 2;

  EulerTools(const double& gmma) // air or water at roughly 20 deg Cels.: gmma = 1.4
      : gamma_(gmma)
  {
  }

  /// \note While the indices may seem silly, using them helps to reduce typo induced errors.
  /// \name To convert between primitive and conservative variables.
  /// \{

  static size_t density_index()
  {
    return 0;
  }

  static std::array<size_t, d> velocity_indices()
  {
    std::array<size_t, d> indices;
    if (d == 1) {
      indices[0] = 1;
    } else if (d == 2) {
      indices[0] = 1;
      indices[1] = 2;
    } else {
      indices[0] = 1;
      indices[1] = 2;
      indices[2] = 3;
    }
    return indices;
  } // ... velocity_indices(...)

  static size_t energy_index()
  {
    return m - 1;
  }

  static size_t pressure_index()
  {
    return m - 1;
  }

  XT::Common::FieldVector<R, 1> density(const FieldVector<R, m>& w) const
  {
    return w[density_index()];
  }

  XT::Common::FieldVector<R, d> velocity(const FieldVector<R, m>& w) const
  {
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = w[velocity_indices()[ii]] / w[density_index()];
    return v;
  }

  XT::Common::FieldVector<R, 1> energy(const FieldVector<R, m>& w) const
  {
    return w[energy_index()];
  }

  XT::Common::FieldVector<R, 1>
  energy(const FieldVector<R, 1>& rho, const FieldVector<R, d>& v, const FieldVector<R, 1>& p) const
  {
    return p / (gamma_ - 1.) + 0.5 * rho * v.two_norm2();
  }

  XT::Common::FieldVector<R, m> primitive(const FieldVector<R, m>& w) const
  {
    XT::Common::FieldVector<R, m> rho_v_p;
    // density
    auto rho = density(w);
    rho_v_p[density_index()] = rho;
    // velocity
    auto v = velocity(w);
    for (size_t ii = 0; ii < d; ++ii)
      rho_v_p[velocity_indices()[ii]] = v[ii];
    // pressure
    rho_v_p[pressure_index()] = (gamma_ - 1.) * (energy(w) - 0.5 * rho * v.two_norm2());
    return rho_v_p;
  } // ... primitive(...)

  std::tuple<XT::Common::FieldVector<R, 1>, XT::Common::FieldVector<R, d>, XT::Common::FieldVector<R, 1>>
  primitives(const FieldVector<R, m>& w) const
  {
    // density
    auto rho = density(w);
    // velocity
    auto v = velocity(w);
    // pressure
    XT::Common::FieldVector<R, 1> p = (gamma_ - 1.) * (energy(w) - 0.5 * rho * v.two_norm2());
    return {std::move(rho), std::move(v), std::move(p)};
  } // ... primitive(...)

  XT::Common::FieldVector<R, 1> pressure(const FieldVector<R, m>& w) const
  {
    return std::get<2>(primitives(w));
  }

  XT::Common::FieldVector<R, m>
  conservative(const FieldVector<R, 1>& rho, const XT::Common::FieldVector<R, d>& v, const FieldVector<R, 1>& p) const
  { // the `v * rho[0]` only works if v is a XT::Common::FieldVector
    return XT::Common::hstack(rho, v * rho[0], energy(rho, v, p));
  }

  /// \}
  /// \name To compute some relevant quantities
  /// \{

  R speed_of_sound2(const FieldVector<R, m>& w) const
  {
    return gamma_ * pressure(w) / density(w);
  }

  R speed_of_sound2(const FieldVector<R, 1>& rho, const FieldVector<R, 1>& p) const
  {
    return gamma_ * p / rho;
  }

  R speed_of_sound(const FieldVector<R, m>& w) const
  {
    return std::sqrt(speed_of_sound2(w));
  }

  R speed_of_sound(const FieldVector<R, 1>& rho, const FieldVector<R, 1>& p) const
  {
    return std::sqrt(speed_of_sound2(rho, p));
  }

  R mach_number(const FieldVector<R, m>& w) const
  {
    const auto rho_v_p = primitives(w);
    const auto& rho = std::get<0>(rho_v_p);
    const auto& v = std::get<1>(rho_v_p);
    const auto& p = std::get<2>(rho_v_p);
    return v.two_norm() / speed_of_sound(rho, p);
  }

  R mach_number(const FieldVector<R, 1>& rho, const FieldVector<R, d>& v, const FieldVector<R, 1>& p) const
  {
    return v.two_norm() / speed_of_sound(rho, p);
  }

  R enthalpy(const FieldVector<R, m>& w) const
  {
    return (energy(w) + pressure(w)) / density(w);
  }

  R enthalpy(const FieldVector<R, 1>& rho, const FieldVector<R, 1>& E, const FieldVector<R, 1>& p) const
  {
    return (E + p) / rho;
  }

  /// \}
  /// \name To access the flux f and its jacobian.
  /// \{

  int flux_order() const
  {
    return 4;
  }

  XT::Common::FieldMatrix<R, d, m> flux(const FieldVector<R, m>& w) const
  {
    // extract
    const auto E = energy(w);
    const auto rho_v_p = primitives(w);
    const auto& rho = std::get<0>(rho_v_p);
    const auto& v = std::get<1>(rho_v_p);
    const auto& p = std::get<2>(rho_v_p);
    // compute
    XT::Common::FieldMatrix<R, d, m> ret;
    for (size_t ss = 0; ss < d; ++ss) {
      auto& f_s = ret[ss];
      f_s[0] = rho * v[ss];
      for (size_t ii = 0; ii < d; ++ii)
        f_s[1 + ii] = rho * v[ii] * v[ss] + (ss == ii ? 1 : 0) * p;
      f_s[m - 1] = (E + p) * v[ss];
    }
    return ret;
  } // ... flux(...)

  /**
   * \brief The Euler flux at impermeable walls (e.g., v*n = 0).
   * \sa    [DF2015, p. 414, (8.58)]
   */
  XT::Common::FieldVector<R, m> flux_at_impermeable_walls(const FieldVector<R, m>& w,
                                                          const XT::Common::FieldVector<double, d>& n) const
  {
    const auto pn = n * pressure(w)[0]; // this * only works if n is a XT::Common::FieldVector
    XT::Common::FieldVector<R, m> ret(0.);
    for (size_t ii = 0; ii < d; ++ii)
      ret[velocity_indices()[ii]] = pn[ii];
    return ret;
  }

  XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d> flux_jacobian(const FieldVector<R, m>& w) const
  {
    // extract
    const auto rho = density(w);
    const auto v = velocity(w);
    const auto E = energy(w);
    // some constants
    const auto gamma_1 = gamma_ - 1.;
    const auto vnorm2 = v.two_norm2();
    const auto ek = 0.5 * vnorm2; // kinetic energy
    // compute
    XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d> ret;
    if (d == 1) {
      auto& jacobian_f_0 = ret[0];
      jacobian_f_0[0] = {0., 1., 0.};
      jacobian_f_0[1] = {gamma_1 * ek - v[0] * v[0], (3. - gamma_) * v[0], gamma_1};
      jacobian_f_0[2] = {v[0] * (gamma_1 * vnorm2 - (gamma_ * E) / rho),
                         ((gamma_ * E) / rho) - gamma_1 * v[0] * v[0] - gamma_1 * ek,
                         gamma_ * v[0]};
    } else if (d == 2) {
      // f_0
      auto& jacobian_f_0 = ret[0];
      jacobian_f_0[0] = {0., 1., 0., 0.};
      jacobian_f_0[1][0] = gamma_1 * ek - v[0] * v[0];
      jacobian_f_0[1][1] = (3. - gamma_) * v[0];
      jacobian_f_0[1][2] = -1. * gamma_1 * v[1];
      jacobian_f_0[1][3] = gamma_1;
      jacobian_f_0[2][0] = -1. * v[0] * v[1];
      jacobian_f_0[2][1] = v[1];
      jacobian_f_0[2][2] = v[0];
      jacobian_f_0[2][3] = 0.;
      jacobian_f_0[3][0] = v[0] * (gamma_1 * vnorm2 - (gamma_ * E) / rho);
      jacobian_f_0[3][1] = ((gamma_ * E) / rho) - gamma_1 * v[0] * v[0] - gamma_1 * ek;
      jacobian_f_0[3][2] = -1. * gamma_1 * v[0] * v[1];
      jacobian_f_0[3][3] = gamma_ * v[0];
      // f_1
      auto& jacobian_f_1 = ret[1];
      jacobian_f_1[0] = {0., 0., 1., 0.};
      jacobian_f_1[1][0] = -1. * v[0] * v[1];
      jacobian_f_1[1][1] = v[1];
      jacobian_f_1[1][2] = v[0];
      jacobian_f_1[1][3] = 0.;
      jacobian_f_1[2][0] = 0.5 * gamma_1 * vnorm2 - v[1] * v[1];
      jacobian_f_1[2][1] = -1. * gamma_1 * v[0];
      jacobian_f_1[2][2] = (3. - gamma_) * v[1];
      jacobian_f_1[2][3] = gamma_1;
      jacobian_f_1[3][0] = v[1] * (gamma_1 * vnorm2 - ((gamma_ * E) / rho));
      jacobian_f_1[3][1] = -1. * gamma_1 * v[0] * v[1];
      jacobian_f_1[3][2] = ((gamma_ * E) / rho) - gamma_1 * v[1] * v[1] - gamma_1 * ek;
      jacobian_f_1[3][3] = gamma_ * v[1];
    } else {
      DUNE_THROW(NotImplemented, "Yet, copy these from DF2015!");
    }
    return ret;
  } // ... flux_jacobian(...)

  /// \}
  /// \{
  /// \name To access the eigendecomposition of the flux jacobian

  /**
   * \sa [Kröner, 1997, p. 387] for 1d and 2d
   */
  XT::Common::FieldVector<R, m> eigenvalues_flux_jacobian(const FieldVector<R, m>& w,
                                                          const FieldVector<double, d>& n) const
  {
    // extract
    const auto rho_v_p = primitives(w);
    const auto& rho = std::get<0>(rho_v_p);
    const auto& v = std::get<1>(rho_v_p);
    const auto& p = std::get<2>(rho_v_p);
    // some constants
    const auto a = speed_of_sound(rho, p);
    const auto vn = v * n;
    // compute
    XT::Common::FieldVector<R, m> ret;
    if (d == 1) {
      // in 1d, we use \lambda_1, \lambda_3, \lambda_4 from the eigendecomposition from [Kröner, 1997, p. 387]
      ret[0] = vn;
      ret[1] = vn + a;
      ret[2] = vn - a;
    } else if (d == 2) {
      // in 2d, we use the eigendecomposition from [Kröner, 1997, p. 387] as is
      ret[0] = vn;
      ret[1] = vn;
      ret[2] = vn + a;
      ret[3] = vn - a;
    } else {
      DUNE_THROW(NotImplemented, "Yet, copy these from Kröner, pp. 391-392!");
    }
    return ret;
  } // ... eigenvalues_flux_jacobian(...)

  /**
   * \sa [Kröner, 1997, p. 387, M T] for 1d and 2d
   */
  XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobian(const FieldVector<R, m>& w,
                                                              const FieldVector<double, d>& n) const
  {
    // extract
    const auto rho_v_p = primitives(w);
    const auto& rho = std::get<0>(rho_v_p);
    const auto& v = std::get<1>(rho_v_p);
    const auto& p = std::get<2>(rho_v_p);
    const auto E = energy(w);
    // some constants
    const auto a = speed_of_sound(rho, p);
    const auto H = enthalpy(rho, E, p);
    const auto rho_over_2a = rho / (2 * a);
    const auto ek = 0.5 * v.two_norm2(); // kinetic energy
    const auto vn = v * n;
    // compute
    XT::Common::FieldMatrix<R, m, m> eigenvectors;
    if (d == 1) {
      // in 1d, we use \lambda_1, \lambda_3, \lambda_4 from the eigendecomposition from [Kröner, 1997, p. 387]
      // this corresponds to striking row 3 and column 2 in MT
      eigenvectors[0][0] = 1.;
      eigenvectors[0][1] = rho_over_2a;
      eigenvectors[0][2] = rho_over_2a;
      eigenvectors[1][0] = v[0];
      eigenvectors[1][1] = rho_over_2a * (v[0] + a * n[0]);
      eigenvectors[1][2] = rho_over_2a * (v[0] - a * n[0]);
      eigenvectors[2][0] = ek;
      eigenvectors[2][1] = rho_over_2a * (H + a * vn);
      eigenvectors[2][2] = rho_over_2a * (H - a * vn);
    } else if (d == 2) {
      // in 2d, we use the eigendecomposition from [Kröner, 1997, p. 387] as is
      // this corresponds to MT
      eigenvectors[0][0] = 1.;
      eigenvectors[0][1] = 0.;
      eigenvectors[0][2] = rho_over_2a;
      eigenvectors[0][3] = rho_over_2a;
      eigenvectors[1][0] = v[0];
      eigenvectors[1][1] = rho * n[1];
      eigenvectors[1][2] = rho_over_2a * (v[0] + a * n[0]);
      eigenvectors[1][3] = rho_over_2a * (v[0] - a * n[0]);
      eigenvectors[2][0] = v[1];
      eigenvectors[2][1] = -rho * n[0];
      eigenvectors[2][2] = rho_over_2a * (v[1] + a * n[1]);
      eigenvectors[2][3] = rho_over_2a * (v[1] - a * n[1]);
      eigenvectors[3][0] = ek;
      eigenvectors[3][1] = rho * (v[0] * n[1] - v[1] * n[0]);
      eigenvectors[3][2] = rho_over_2a * (H + a * vn);
      eigenvectors[3][3] = rho_over_2a * (H - a * vn);
    } else {
      DUNE_THROW(NotImplemented, "Yet, copy these from Kröner, pp. 391-392!");
    }
    return eigenvectors;
  } // ... eigenvectors_flux_jacobian(...)

  /**
   * \sa [Kröner, 1997, p. 387, (M T)^{-1}] (with q = 1 in the top left entry) for 1d and 2d
   */
  XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobian(const FieldVector<R, m>& w,
                                                                  const FieldVector<double, d>& n) const
  {
    // extract
    const auto rho_v_p = primitives(w);
    const auto& rho = std::get<0>(rho_v_p);
    const auto& v = std::get<1>(rho_v_p);
    const auto& p = std::get<2>(rho_v_p);
    // some constants
    const auto a = speed_of_sound(rho, p);
    const auto M = mach_number(rho, v, p);
    const auto vn = v * n;
    const auto gamma_1 = gamma_ - 1.;
    // compute
    XT::Common::FieldMatrix<R, m, m> eigenvectors_inv;
    if (d == 1) {
      // in 1d, we use \lambda_1, \lambda_3, \lambda_4 from the eigendecomposition from [Kröner, 1997, p. 387]
      // this corresponds to striking row 3 and column 2 in (MT)^{-1}
      eigenvectors_inv[0][0] = 1. - (gamma_1 / 2.) * M * M;
      eigenvectors_inv[0][1] = gamma_1 * v[0] / (a * a);
      eigenvectors_inv[0][2] = -gamma_1 / (a * a);
      eigenvectors_inv[1][0] = (a / rho) * ((gamma_1 / 2.) * M * M - vn / a);
      eigenvectors_inv[1][1] = (1. / rho) * (n[0] - gamma_1 * (v[0] / a));
      eigenvectors_inv[1][2] = gamma_1 / (rho * a);
      eigenvectors_inv[2][0] = (a / rho) * ((gamma_1 / 2.) * M * M + vn / a);
      eigenvectors_inv[2][1] = (-1. / rho) * (n[0] + gamma_1 * (v[0] / a));
      eigenvectors_inv[2][2] = gamma_1 / (rho * a);
    } else if (d == 2) {
      // in 2d, we use the eigendecomposition from [Kröner, 1997, p. 387] as is
      // this corresponds to (MT)^{-1}
      eigenvectors_inv[0][0] = 1. - (gamma_1 / 2.) * M * M;
      eigenvectors_inv[0][1] = gamma_1 * v[0] / (a * a);
      eigenvectors_inv[0][2] = gamma_1 * v[1] / (a * a);
      eigenvectors_inv[0][3] = -gamma_1 / (a * a);
      eigenvectors_inv[1][0] = (1. / rho) * (v[1] * n[0] - v[0] * n[1]);
      eigenvectors_inv[1][1] = n[1] / rho;
      eigenvectors_inv[1][2] = -n[0] / rho;
      eigenvectors_inv[1][3] = 0.;
      eigenvectors_inv[2][0] = (a / rho) * ((gamma_1 / 2.) * M * M - vn / a);
      eigenvectors_inv[2][1] = (1. / rho) * (n[0] - gamma_1 * (v[0] / a));
      eigenvectors_inv[2][2] = (1. / rho) * (n[1] - gamma_1 * (v[1] / a));
      eigenvectors_inv[2][3] = gamma_1 / (rho * a);
      eigenvectors_inv[3][0] = (a / rho) * ((gamma_1 / 2.) * M * M + vn / a);
      eigenvectors_inv[3][1] = (-1. / rho) * (n[0] + gamma_1 * (v[0] / a));
      eigenvectors_inv[3][2] = (-1. / rho) * (n[1] + gamma_1 * (v[1] / a));
      eigenvectors_inv[3][3] = gamma_1 / (rho * a);
    } else {
      DUNE_THROW(NotImplemented, "Yet, copy these from Kröner, pp. 391-392!");
    }
    return eigenvectors_inv;
  } // ... eigenvectors_inv_flux_jacobian(...)

  /// \}

  template <class E, class GL>
  std::enable_if_t<XT::Grid::is_layer<GL>::value, void>
  visualize(const XT::Functions::GridFunctionInterface<E, m, 1, R>& u_conservative,
            const GL& grid_layer,
            const std::string filename_prefix = "",
            const std::string filename_suffix = "",
            const bool subsampling = false) const
  {
    const std::string prefix = filename_prefix.empty() ? "" : filename_prefix + "_";
    const std::string suffix = filename_suffix.empty() ? "" : "_" + filename_suffix;
    XT::Functions::make_sliced_function<1>(u_conservative, {density_index()}, "density")
        .visualize(grid_layer, prefix + "density" + suffix, subsampling);
    XT::Functions::make_sliced_function<d>(u_conservative, velocity_indices(), "density_times_velocity")
        .visualize(grid_layer, prefix + "density_times_velocity" + suffix, subsampling);
    XT::Functions::make_sliced_function<1>(u_conservative, {energy_index()}, "energy")
        .visualize(grid_layer, prefix + "energy" + suffix, subsampling);
    const auto u_primitive =
        XT::Functions::make_transformed_function<m, 1, R>(u_conservative, [&](const auto& w) { return primitive(w); });
    XT::Functions::make_sliced_function<d>(u_primitive, velocity_indices(), "velocity")
        .visualize(grid_layer, prefix + "velocity" + suffix, subsampling);
    XT::Functions::make_sliced_function<1>(u_primitive, {pressure_index()}, "pressure")
        .visualize(grid_layer, prefix + "pressure" + suffix, subsampling);
  } // ... visualize(...)

private:
  const double gamma_;
}; // class EulerTools


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_EULER_HH
