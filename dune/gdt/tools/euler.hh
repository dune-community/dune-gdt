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
 */
template <size_t d, class R = double>
class EulerTools
{
public:
  static const constexpr size_t m = d + 2;

  EulerTools(const double& gmma) // air or water at roughly 20 deg Cels.: gmma = 1.4
      : gamma_(gmma)
  {
  }

private:
  template <class VectorType>
  XT::Common::FieldVector<R, m> convert_to_primitive(const VectorType& conservative_variables) const
  {
    // extract
    const auto& rho = conservative_variables[0];
    FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[ii + 1] / rho;
    const auto& e = conservative_variables[m - 1];
    // convert
    XT::Common::FieldVector<R, m> primitive_variables;
    // * density
    primitive_variables[0] = rho;
    // * velocity
    for (size_t ii = 0; ii < d; ++ii)
      primitive_variables[ii + 1] = v[ii];
    // * pressure
    primitive_variables[m - 1] = (gamma_ - 1.) * (e - 0.5 * rho * v.two_norm2());
    return primitive_variables;
  }

public:
  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11)]
   */
  XT::Common::FieldVector<R, m> to_primitive(const FieldVector<R, m>& conservative_variables) const
  {
    return convert_to_primitive(conservative_variables);
  }

  template <class V>
  XT::Common::FieldVector<R, m> to_primitive(const XT::LA::VectorInterface<V>& conservative_variables) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return convert_to_primitive(conservative_variables);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, m> to_conservative(const FieldVector<R, m>& primitive_variables) const
  {
    // extract
    const auto& rho = primitive_variables[0];
    FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = primitive_variables[ii + 1];
    const auto& p = primitive_variables[m - 1];
    // convert
    XT::Common::FieldVector<R, m> conservative_variables;
    // * density
    conservative_variables[0] = rho;
    // * density times velocity component
    for (size_t ii = 0; ii < d; ++ii)
      conservative_variables[1 + ii] = rho * v[ii];
    // * energy
    conservative_variables[m - 1] = p / (gamma_ - 1.) + 0.5 * rho * v.two_norm2();
    return conservative_variables;
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, 1> density_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    return conservative_variables[0];
  }

  template <class V>
  XT::Common::FieldVector<R, 1>
  density_from_conservative(const XT::LA::VectorInterface<V>& conservative_variables) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return conservative_variables[0];
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, 1> density_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    return primitive_variables[0];
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, d>
  density_times_velocity_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[1 + ii];
    return v;
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, d>
  density_times_velocity_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    const auto& rho = primitive_variables[0];
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = rho * primitive_variables[1 + ii];
    return v;
  }

private:
  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  template <class VectorType>
  XT::Common::FieldVector<R, d> extract_velocity_from_conservative(const VectorType& conservative_variables) const
  {
    const auto& rho = conservative_variables[0];
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[1 + ii] / rho;
    return v;
  }

public:
  XT::Common::FieldVector<R, d> velocity_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    return extract_velocity_from_conservative(conservative_variables);
  }

  template <class V>
  XT::Common::FieldVector<R, d>
  velocity_from_conservative(const XT::LA::VectorInterface<V>& conservative_variables) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return extract_velocity_from_conservative(conservative_variables);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, d> velocity_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = primitive_variables[1 + ii];
    return v;
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, 1> energy_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    return conservative_variables[m - 1];
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, 1> energy_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    const auto& rho = primitive_variables[0];
    FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = primitive_variables[ii + 1];
    const auto& p = primitive_variables[m - 1];
    return p / (gamma_ - 1.) + 0.5 * rho * v.two_norm2();
  }

private:
  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  template <class VectorType>
  XT::Common::FieldVector<R, 1> extract_pressure_from_conservative(const VectorType& conservative_variables) const
  {
    const auto& rho = conservative_variables[0];
    FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[ii + 1] / rho;
    const auto& e = conservative_variables[m - 1];
    return (gamma_ - 1.) * (e - 0.5 * rho * v.two_norm2());
  }

public:
  XT::Common::FieldVector<R, 1> pressure_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    return extract_pressure_from_conservative(conservative_variables);
  }

  template <class V>
  XT::Common::FieldVector<R, 1>
  pressure_from_conservative(const XT::LA::VectorInterface<V>& conservative_variables) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return extract_pressure_from_conservative(conservative_variables);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 404, (8.11) or p. 421, (8.89)]
   */
  XT::Common::FieldVector<R, 1> pressure_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    return primitive_variables[m - 1];
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 403, (8.10)]
   */
  XT::Common::FieldMatrix<R, d, m> flux(const FieldVector<R, m>& conservative_variables) const
  {
    const auto primitive_variables = to_primitive(conservative_variables);
    const auto& rho = conservative_variables[0];
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v = primitive_variables[ii + 1];
    const auto& e = conservative_variables[m - 1];
    const auto& p = primitive_variables[m - 1];
    XT::Common::FieldMatrix<R, d, m> ret;
    for (size_t ss = 0; ss < d; ++ss) {
      auto& f_s = ret[ss];
      f_s[0] = rho * v[ss];
      for (size_t ii = 0; ii < d; ++ii)
        f_s[1 + ii] = rho * v[ii] * v[ss] + (ss == ii ? 1 : 0) * p;
      f_s[m - 1] = (e + p) * v[ss];
    }
    return ret;
  } // ... flux(...)

  int flux_order() const
  {
    return 4;
  }

private:
  /**
   * \brief The Euler flux at impermeable walls (e.g., v*n = 0).
   * \sa    [DF2015, p. 414, (8.58)]
   */
  template <class VectorType, class D>
  XT::Common::FieldVector<R, m> compute_flux_at_impermeable_walls(const VectorType& conservative_variables,
                                                                  const FieldVector<D, d>& normal) const
  {
    const auto pressure = to_primitive(conservative_variables)[m - 1];
    const auto tmp = normal * pressure;
    XT::Common::FieldVector<R, m> ret(0.);
    for (size_t ii = 0; ii < d; ++ii)
      ret[ii + 1] = tmp[ii];
    return ret;
  } // ... flux_at_impermeable_walls(...)

public:
  template <class D>
  XT::Common::FieldVector<R, m> flux_at_impermeable_walls(const FieldVector<R, m>& conservative_variables,
                                                          const FieldVector<D, d>& normal) const
  {
    return compute_flux_at_impermeable_walls(conservative_variables, normal);
  }

  template <class V, class D>
  XT::Common::FieldVector<R, m> flux_at_impermeable_walls(const XT::LA::VectorInterface<V>& conservative_variables,
                                                          const FieldVector<D, d>& normal) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return compute_flux_at_impermeable_walls(conservative_variables, normal);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 405, (8.18 - 8.19)] for the 2d case
   */
  XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
  flux_jacobian(const FieldVector<R, m>& conservative_variables) const
  {
    return dim_switch<>::flux_jacobian(gamma_, conservative_variables);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 405, (8.20)] for the 2d case
   */
  XT::Common::FieldMatrix<R, m, m> flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                      const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::flux_jacobi_matrix(gamma_, conservative_variables, normal);
  }

  /**
   * \sa [Kröner, 1997, p. 387, \lambda_1, \dots, \lambda_4] for the 2d case
   */
  XT::Common::FieldVector<R, m> eigenvalues_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                               const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvalues_flux_jacobi_matrix(gamma_, conservative_variables, normal);
  }

  XT::Common::FieldMatrix<R, m, m> eigenvaluematrix_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                       const FieldVector<double, d>& normal) const
  {
    const auto evs = eigenvalues_flux_jacobi_matrix(conservative_variables, normal);
    XT::Common::FieldMatrix<R, m, m> ret(0.);
    for (size_t ii = 0; ii < m; ++ii)
      ret[ii][ii] = evs[ii];
    return ret;
  }

  /**
   * \sa [Kröner, 1997, p. 387, M T] for the 2d case
   */
  XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                   const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvectors_flux_jacobi_matrix(gamma_, conservative_variables, normal);
  }

  /**
   * \sa [Kröner, 1997, p. 387, (M T)^{-1}] for the 2d case (with q = 1 in the top left entry)
   */
  XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                       const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvectors_inv_flux_jacobi_matrix(gamma_, conservative_variables, normal);
  }

private:
  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 403, (8.7)] for the 2d case
   */
  template <class VectorType>
  R compute_speed_of_sound_from_conservative(const VectorType& conservative_variables) const
  {
    const auto rho = density_from_conservative(conservative_variables);
    const auto p = pressure_from_conservative(conservative_variables);
    return std::sqrt((gamma_ * p) / rho);
  }

public:
  R speed_of_sound_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    return compute_speed_of_sound_from_conservative(conservative_variables);
  }

  template <class V>
  R speed_of_sound_from_conservative(const XT::LA::VectorInterface<V>& conservative_variables) const
  {
    DUNE_THROW_IF(conservative_variables.size() != m,
                  XT::Common::Exceptions::shapes_do_not_match,
                  "conservative_variables.size() = " << conservative_variables.size() << "\n   m = " << m);
    return compute_speed_of_sound_from_conservative(conservative_variables);
  }

  /**
   * \sa [Dolejsi, Feistauer, 2016, p. 403, (8.7)] for the 2d case
   */
  R speed_of_sound_from_primitive(const FieldVector<R, m>& primitive_variables) const
  {
    const auto rho = density_from_primiitve(primitive_variables);
    const auto p = pressure_from_primitive(primitive_variables);
    return std::sqrt((gamma_ * p) / rho);
  }

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
    XT::Functions::make_sliced_function<1>(u_conservative, {0}, "density")
        .visualize(grid_layer, prefix + "density" + suffix, subsampling);
    XT::Functions::make_sliced_function<d>(u_conservative, dim_switch<>::velocity_indices(), "density_times_velocity")
        .visualize(grid_layer, prefix + "density_times_velocity" + suffix, subsampling);
    XT::Functions::make_sliced_function<1>(u_conservative, {m - 1}, "energy")
        .visualize(grid_layer, prefix + "energy" + suffix, subsampling);
    const auto u_primitive = XT::Functions::make_transformed_function<m, 1, R>(
        u_conservative, [&](const auto& cons) { return to_primitive(cons); });
    XT::Functions::make_sliced_function<d>(u_primitive, dim_switch<>::velocity_indices(), "velocity")
        .visualize(grid_layer, prefix + "velocity" + suffix, subsampling);
    XT::Functions::make_sliced_function<1>(u_primitive, {m - 1}, "pressure")
        .visualize(grid_layer, prefix + "pressure" + suffix, subsampling);
  } // ... visualize(...)

private:
  template <size_t d_ = d, typename anything = void>
  struct dim_switch;

  template <typename anything>
  struct dim_switch<1, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1};
    }

    static XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
    flux_jacobian(const double& gamma, const FieldVector<R, m>& conservative_variables)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const auto& v = primitive_variables[1];
      const auto& e = conservative_variables[2];
      XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d> ret;
      auto& jacobian_f = ret[0];
      jacobian_f[0][0] = 0.;
      jacobian_f[0][1] = 1.;
      jacobian_f[0][2] = 0.;
      jacobian_f[1][0] = 0.5 * (gamma - 3.) * v * v;
      jacobian_f[1][1] = (3. - gamma) * v;
      jacobian_f[1][2] = gamma - 1.;
      jacobian_f[2][0] = v * ((gamma - 1.) * v * v - gamma * (e / rho));
      jacobian_f[2][1] = gamma * (e / rho) - (3. * (gamma - 1.) / 2.) * v * v;
      jacobian_f[2][2] = gamma * v;
      return ret;
    } // ... flux_jacobian(...)
  }; // struct dim_switch<1, ...>

  template <typename anything>
  struct dim_switch<2, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1, 2};
    }

    /**
     * \sa [Dolejsi, Feistauer, 2016, p. 405, (8.18 - 8.19)]
     */
    static XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
    flux_jacobian(const double& gamma, const FieldVector<R, m>& conservative_variables)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const XT::Common::FieldVector<R, d> v = {primitive_variables[1], primitive_variables[2]};
      const auto& e = conservative_variables[m - 1];
      XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d> ret;
      const auto gamma_1 = gamma - 1.;
      // f_0
      auto& jacobian_f_0 = ret[0];
      jacobian_f_0[0] = {0., 1., 0., 0.};
      jacobian_f_0[1][0] = 0.5 * gamma_1 * v.two_norm2() - v[0] * v[0];
      jacobian_f_0[1][1] = (3. - gamma) * v[0];
      jacobian_f_0[1][2] = -1. * gamma_1 * v[1];
      jacobian_f_0[1][3] = gamma_1;
      jacobian_f_0[2][0] = -1. * v[0] * v[1];
      jacobian_f_0[2][1] = v[1];
      jacobian_f_0[2][2] = v[0];
      jacobian_f_0[2][3] = 0.;
      jacobian_f_0[3][0] = v[0] * (gamma_1 * v.two_norm2() - (gamma * e) / rho);
      jacobian_f_0[3][1] = ((gamma * e) / rho) - gamma_1 * v[0] * v[0] - 0.5 * gamma_1 * v.two_norm2();
      jacobian_f_0[3][2] = -1. * gamma_1 * v[0] * v[1];
      jacobian_f_0[3][3] = gamma * v[0];
      // f_1
      auto& jacobian_f_1 = ret[1];
      jacobian_f_1[0] = {0., 0., 1., 0.};
      jacobian_f_1[1][0] = -1. * v[0] * v[1];
      jacobian_f_1[1][1] = v[1];
      jacobian_f_1[1][2] = v[0];
      jacobian_f_1[1][3] = 0.;
      jacobian_f_1[2][0] = 0.5 * gamma_1 * v.two_norm2() - v[1] * v[1];
      jacobian_f_1[2][1] = -1. * gamma_1 * v[0];
      jacobian_f_1[2][2] = (3. - gamma) * v[1];
      jacobian_f_1[2][3] = gamma_1;
      jacobian_f_1[3][0] = v[1] * (gamma_1 * v.two_norm2() - ((gamma * e) / rho));
      jacobian_f_1[3][1] = -1. * gamma_1 * v[0] * v[1];
      jacobian_f_1[3][2] = ((gamma * e) / rho) - gamma_1 * v[1] * v[1] - 0.5 * gamma_1 * v.two_norm2();
      jacobian_f_1[3][3] = gamma * v[1];
      return ret;
    } // ... flux_jacobian(...)

    /**
     * \sa [Dolejsi, Feistauer, 2016, p. 405, (8.20)]
     */
    static XT::Common::FieldMatrix<R, m, m> flux_jacobi_matrix(const double& gamma,
                                                               const FieldVector<R, m>& conservative_variables,
                                                               const FieldVector<double, d>& normal)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const XT::Common::FieldVector<R, d> v = {primitive_variables[1], primitive_variables[2]};
      const auto& e = conservative_variables[3];
      const auto gamma_1 = gamma - 1;
      const auto gamma_2 = gamma - 2;
      const auto G = gamma * e / rho - 0.5 * gamma_1 * v.two_norm2();
      XT::Common::FieldMatrix<R, m, m> P(0.);
      P[0][0] = 0;
      P[0][1] = normal[0];
      P[0][2] = normal[1];
      P[0][3] = 0;

      P[1][0] = 0.5 * gamma_1 * v.two_norm2() * normal[0] - v[0] * (v * normal);
      P[1][1] = -gamma_2 * v[0] * normal[0] + v * normal;
      P[1][2] = v[0] * normal[1] - gamma_1 * v[1] * normal[0];
      P[1][3] = gamma_1 * normal[0];

      P[2][0] = 0.5 * gamma_1 * v.two_norm2() * normal[1] - v[1] * (v * normal);
      P[2][1] = v[1] * normal[0] - gamma_1 * v[0] * normal[1];
      P[2][2] = -gamma_2 * v[1] * normal[1] + v * normal;
      P[2][3] = gamma_1 * normal[1];

      P[3][0] = (gamma_1 * v.two_norm2() - gamma * e / rho) * (v * normal);
      P[3][1] = G * normal[0] - gamma_1 * v[0] * (v * normal);
      P[3][2] = G * normal[1] - gamma_1 * v[1] * (v * normal);
      P[3][3] = gamma * (v * normal);
      return P;
    } // ... flux_jacobi_matrix(...)

    /**
     * \sa [Kröner, 1997, p. 387, \lambda_1, \dots, \lambda_4]
     */
    static XT::Common::FieldVector<R, m> eigenvalues_flux_jacobi_matrix(const double& gamma,
                                                                        const FieldVector<R, m>& conservative_variables,
                                                                        const FieldVector<double, d>& normal)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const auto& p = primitive_variables[3];
      const auto a = std::sqrt(gamma * p / rho);
      const auto v_times_n = XT::Common::FieldVector<R, d>({primitive_variables[1], primitive_variables[2]}) * normal;
      return {v_times_n, v_times_n, v_times_n + a, v_times_n - a};
    } // ... eigenvalues_flux_jacobi_matrix(...)

    /**
     * \sa [Kröner, 1997, p. 387, M T]
     */
    static XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobi_matrix(
        const double& gamma, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& normal)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const auto& v_0 = primitive_variables[1];
      const auto& v_1 = primitive_variables[2];
      const auto v_abs_2 = v_0 * v_0 + v_1 * v_1;
      const auto& e = conservative_variables[3];
      const auto& p = primitive_variables[3];
      const auto a = std::sqrt(gamma * p / rho);
      const auto H = (e + p) / rho;
      const auto rho_over_2a = rho / (2 * a);
      const auto v_times_n = v_0 * normal[0] + v_1 * normal[1];
      XT::Common::FieldMatrix<R, m, m> eigenvectors(0.);
      eigenvectors[0] = {1., 0., rho_over_2a, rho_over_2a};
      eigenvectors[1] = {
          v_0, rho * normal[1], rho_over_2a * (v_0 + a * normal[0]), rho_over_2a * (v_0 - a * normal[0])};
      eigenvectors[2] = {
          v_1, -rho * normal[0], rho_over_2a * (v_1 + a * normal[1]), rho_over_2a * (v_1 - a * normal[1])};
      eigenvectors[3] = {v_abs_2 / 2.,
                         rho * (v_0 * normal[1] - v_1 * normal[0]),
                         rho_over_2a * (H + a * v_times_n),
                         rho_over_2a * (H - a * v_times_n)};
      return eigenvectors;
    } // ... eigenvectors_flux_jacobi_matrix(...)

    /**
     * \sa [Kröner, 1997, p. 387, (M T)^{-1}] (with q = 1 in the top left entry)
     */
    static XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobi_matrix(
        const double& gamma, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& normal)
    {
      const auto primitive_variables = EulerTools<d, R>(gamma).to_primitive(conservative_variables);
      const auto& rho = conservative_variables[0];
      const auto& v_0 = primitive_variables[1];
      const auto& v_1 = primitive_variables[2];
      const auto v_abs_2 = v_0 * v_0 + v_1 * v_1;
      const auto v_abs = std::sqrt(v_abs_2);
      const auto& p = primitive_variables[3];
      const auto a = std::sqrt(gamma * p / rho);
      const auto M = v_abs / a;
      const auto v_times_n = v_0 * normal[0] + v_1 * normal[1];

      XT::Common::FieldMatrix<R, m, m> eigenvectors_inv(0.);
      eigenvectors_inv[0] = {1. - ((gamma - 1.) / 2.) * M * M,
                             (gamma - 1.) * v_0 / (a * a),
                             (gamma - 1.) * v_1 / (a * a),
                             -(gamma - 1.) / (a * a)};
      eigenvectors_inv[1] = {(1. / rho) * (v_1 * normal[0] - v_0 * normal[1]), normal[1] / rho, -normal[0] / rho, 0.};
      eigenvectors_inv[2] = {(a / rho) * (((gamma - 1.) / 2) * M * M - v_times_n / a),
                             (1. / rho) * (normal[0] - (gamma - 1.) * (v_0 / a)),
                             (1. / rho) * (normal[1] - (gamma - 1.) * (v_1 / a)),
                             (gamma - 1.) / (rho * a)};
      eigenvectors_inv[3] = {(a / rho) * (((gamma - 1.) / 2) * M * M + v_times_n / a),
                             (-1. / rho) * (normal[0] + (gamma - 1.) * (v_0 / a)),
                             (-1. / rho) * (normal[1] + (gamma - 1.) * (v_1 / a)),
                             (gamma - 1.) / (rho * a)};
      return eigenvectors_inv;
    } // ... eigenvectors_inv_flux_jacobi_matrix(...)
  }; //  struct dim_switch<2, ...>

  template <typename anything>
  struct dim_switch<3, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1, 2, 3};
    }
  };

  const double gamma_;
}; // class EulerTools


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_EULER_HH
