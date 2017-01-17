// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
#define DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH

#include <cmath>
#include <algorithm>
#include <memory>
#include <unordered_map>

#include <dune/gdt/local/fluxes/interfaces.hh>

namespace Dune {
namespace GDT {


template <class R, int r>
FieldMatrix<R, r, r> unit_matrix()
{
  FieldMatrix<R, r, r> ret(0);
  for (size_t ii = 0; ii < r; ++ii)
    ret[ii][ii] = 1;
  return ret;
}

/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of \mathbf{f}.
 */
template <class GridType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC, size_t num_quad_points>
class EntropyBasedFlux : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedFlux<GridType, E, D, d, R, rangeDim, rC, num_quad_points> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using BasisValuesMatrixType = Dune::FieldMatrix<RangeFieldType, num_quad_points, dimRange>;
  using typename BaseType::FluxRangeType;
  using typename BaseType::FluxJacobianRangeType;

  explicit EntropyBasedFlux(
      const GridType& grid,
      const Dune::QuadratureRule<DomainFieldType, dimDomain>& quadrature,
      const BasisValuesMatrixType& M,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = unit_matrix<RangeFieldType, dimRange>(),
      const std::string name = static_id())
    : id_set_(grid.localIdSet())
    , quadrature_(quadrature)
    , M_(M)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , T_minus_one_(T_minus_one)
    , name_(name)
  {
    assert(quadrature.size() == num_quad_points);
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& /*u*/, const E& /*entity*/, const DomainType& /*x_local*/, const double /*t*/ = 1) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  virtual FluxRangeType
  evaluate(const RangeType& u, const E& entity, const DomainType& /*x_local*/, const double t = 0) const
  {
    // get initial multiplier and basis matrix from last time step
    RangeType alpha_iso(0);
    alpha_iso[0] = std::log(u[0] / 2.);
    const auto id = id_set_.id(entity);
    const auto beta_it = beta_cache_.find(id);
    const auto T_it = T_cache_.find(id);
    BasisValuesMatrixType P_k;

    // define further variables
    bool chol_flag = false;
    RangeType g_k, beta_out, alpha;
    MatrixType T_k;

    const auto r_max = r_sequence_.back();
    for (const auto& r : r_sequence_) {
      RangeType beta_in = beta_it != beta_cache_.end() && t > -99 ? beta_it->second : alpha_iso;
      T_k = T_it != T_cache_.end() && t > -99 ? T_it->second : T_minus_one_;
      // normalize u
      RangeType r_times_u_iso(0);
      r_times_u_iso[0] = r * u[0];
      RangeType v = u;
      v *= 1 - r;
      v += r_times_u_iso;
      //      std::cout << "v ist : " << XT::Common::to_string(v) << std::endl;
      // calculate T_k u
      RangeType v_k;
      T_k.solve(v_k, v);
      // calculate values of basis p = T_k m
      for (size_t ii = 0; ii < num_quad_points; ++ii)
        T_k.solve(P_k[ii], M_[ii]);
      // calculate f_0
      RangeFieldType f_k(0);
      for (size_t ll = 0; ll < quadrature_.size(); ++ll)
        f_k += quadrature_[ll].weight() * std::exp(beta_in * P_k[ll]);
      f_k -= beta_in * v_k;

      for (size_t kk = 0; kk < k_max_; ++kk) {
        //        std::cout << "k: " << kk << " " << ", r: " << r << " " <<  std::endl;
        change_basis(chol_flag, beta_in, v_k, P_k, T_k, g_k, beta_out);
        MatrixType T_k_transp(0);
        for (size_t ii = 0; ii < dimRange; ++ii)
          for (size_t jj = 0; jj < dimRange; ++jj)
            T_k_transp[ii][jj] = T_k[jj][ii];
        if (chol_flag == false && r == r_max)
          DUNE_THROW(Dune::NotImplemented, "Failure to converge!");
        // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
        if ((kk > k_0_ || chol_flag == false) && r < r_max)
          break;
        // calculate current error
        RangeType error(0);
        for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
          auto m = M_[ll];
          RangeType Tinv_m(0);
          T_k.solve(Tinv_m, m);
          m *= std::exp(beta_out * Tinv_m) * quadrature_[ll].weight();
          error += m;
        }
        error -= v;
        // calculate descent direction d_k;
        RangeType d_k = g_k;
        d_k *= -1;
        RangeType T_k_inv_transp_d_k;
        T_k_transp.solve(T_k_inv_transp_d_k, d_k);
        if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.one_norm()) < 1 + epsilon_gamma_) {
          T_k_transp.solve(alpha, beta_out);
          goto outside_all_loops;
        } else {
          RangeFieldType zeta_k = 1;
          beta_in = beta_out;
          // backtracking line search
          while (zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm()) {
            RangeFieldType f(0);
            auto beta_new = d_k;
            beta_new *= zeta_k;
            beta_new += beta_out;
            for (size_t ll = 0; ll < quadrature_.size(); ++ll)
              f += quadrature_[ll].weight() * std::exp(beta_new * P_k[ll]);
            f -= beta_new * v_k;
            if (XT::Common::FloatCmp::le(f, f_k + xi_ * zeta_k * (g_k * d_k))) {
              beta_in = beta_new;
              f_k = f;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // r loop (Regularization parameter)

    DUNE_THROW(NotImplemented, "Failed to converge");

  outside_all_loops:
    // store values as initial conditions for next time step on this entity
    beta_cache_[id] = beta_out;
    T_cache_[id] = T_k;

    // calculate < \mu m G_\alpha(u) >
    RangeType ret(0);
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& mu = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      auto m = M_[ll];
      m *= mu * std::exp(alpha * m) * weight;
      ret += m;
    }
    return ret;
  } // FluxRangeType evaluate(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  void change_basis(bool& chol_flag,
                    const RangeType& beta_in,
                    RangeType& v_k,
                    BasisValuesMatrixType& P_k,
                    MatrixType& T_k,
                    RangeType& g_k,
                    RangeType& beta_out) const
  {
    MatrixType H(0), tmp(0), L(0);
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& p = P_k[ll];
      // calculate p p^T
      for (size_t ii = 0; ii < dimRange; ++ii) {
        tmp[ii] = p;
        tmp[ii] *= p[ii];
      }
      // multiply p p^T by factor
      tmp *= std::exp(beta_in * p) * quadrature_[ll].weight();
      H += tmp;
    } // quadrature points for loop
    chol_flag = cholesky_L(H, L);
    if (chol_flag == false)
      return;
    const auto P_k_copy = P_k;
    const auto v_k_copy = v_k;
    for (size_t ll = 0; ll < num_quad_points; ++ll)
      L.solve(P_k[ll], P_k_copy[ll]);
    T_k.rightmultiply(L);
    L.mtv(beta_in, beta_out);
    L.solve(v_k, v_k_copy);
    g_k = v_k;
    g_k *= -1;
    // assumes that the first basis function is constant
    for (size_t ll = 0; ll < quadrature_.size(); ++ll)
      g_k[0] += std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight() * P_k[0][0];
  } // void change_basis(...)

  // copied and adapted from dune/geometry/affinegeometry.hh
  template <int size>
  static bool cholesky_L(const FieldMatrix<RangeFieldType, size, size>& H, FieldMatrix<RangeFieldType, size, size>& L)
  {
    for (int ii = 0; ii < size; ++ii) {
      RangeFieldType& rii = L[ii][ii];

      RangeFieldType xDiag = H[ii][ii];
      for (int jj = 0; jj < ii; ++jj)
        xDiag -= L[ii][jj] * L[ii][jj];

      //      std::cout << XT::Common::to_string(H) << std::endl;
      if (XT::Common::FloatCmp::le(xDiag, RangeFieldType(0)))
        return false;

      rii = std::sqrt(xDiag);

      RangeFieldType invrii = RangeFieldType(1) / rii;
      for (int kk = ii + 1; kk < size; ++kk) {
        RangeFieldType x = H[kk][ii];
        for (int jj = 0; jj < ii; ++jj)
          x -= L[ii][jj] * L[kk][jj];
        L[kk][ii] = invrii * x;
      }
    }
    return true;
  }

  const typename GridType::LocalIdSet& id_set_;
  const Dune::QuadratureRule<DomainFieldType, dimDomain>& quadrature_;
  const Dune::FieldMatrix<RangeFieldType, num_quad_points, dimRange> M_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const MatrixType T_minus_one_;
  const std::string name_;
  mutable std::unordered_map<typename GridType::LocalIdSet::IdType, RangeType> beta_cache_;
  mutable std::unordered_map<typename GridType::LocalIdSet::IdType, MatrixType> T_cache_;
};


/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of \mathbf{f}.
 */
template <class GridType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class EntropyBasedFluxHatFunctions : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedFluxHatFunctions<GridType, E, D, d, R, rangeDim, rC> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using typename BaseType::FluxRangeType;
  using typename BaseType::FluxJacobianRangeType;

  explicit EntropyBasedFluxHatFunctions(
      const GridType& grid,
      const RangeType v_points,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const RangeFieldType taylor_tol = 1e-10,
      const std::string name = static_id())
    : id_set_(grid.localIdSet())
    , v_points_(v_points)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , taylor_tol_(taylor_tol)
    , name_(name)
  {
    assert(quadrature.size() == num_quad_points);
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& /*u*/, const E& /*entity*/, const DomainType& /*x_local*/, const double /*t*/ = 1) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  // integrals can be evaluated exactly for hatfunctions
  FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& /*x_local*/, const double t = 0) const
  {
    // get initial multiplier and basis matrix from last time step
    RangeType alpha_iso(1);
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 2.;
    alpha_iso *= std::log(psi_iso);
    const auto id = id_set_.id(entity);
    const auto alpha_it = alpha_cache_.find(id);

    // define further variables
    RangeType g_k, alpha;
    MatrixType H_k;

    // calculate moment vector for isotropic distribution
    RangeType u_iso(0);
    u_iso[0] = v_points_[1] - v_points_[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      u_iso[ii] = v_points_[ii + 1] - v_points_[ii - 1];
    u_iso[dimRange - 1] = v_points_[dimRange - 1] - v_points_[dimRange - 2];
    u_iso *= psi_iso / 2.;

    const auto r_max = r_sequence_.back();
    for (const auto& r : r_sequence_) {
      // get initial alpha
      RangeType alpha_k = alpha_it != alpha_cache_.end() && t > -99 ? alpha_it->second : alpha_iso;
      // normalize u
      RangeType r_times_u_iso(u_iso);
      r_times_u_iso *= r;
      RangeType v = u;
      v *= 1 - r;
      v += r_times_u_iso;

      // calculate f_0
      RangeFieldType f_k(0);
      for (size_t ii = 0; ii < dimRange - 1; ++ii) {
        if (XT::Common::FloatCmp::ne(alpha_k[ii + 1], alpha_k[ii], taylor_tol_))
          f_k += (v_points_[ii + 1] - v_points_[ii]) / (alpha_k[ii + 1] - alpha_k[ii])
                 * (std::exp(alpha_k[ii + 1]) - std::exp(alpha_k[ii]));
        else
          f_k += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_k[ii]);
      }
      f_k -= alpha_k * v;

      for (size_t kk = 0; kk < k_max_; ++kk) {
        // calculate gradient g
        g_k *= 0;
        for (size_t nn = 0; nn < dimRange; ++nn) {
          if (nn > 0) {
            if (XT::Common::FloatCmp::ne(alpha_k[nn], alpha_k[nn - 1], taylor_tol_)) {
              g_k[nn] +=
                  -(v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
                      * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                  + (v_points_[nn] - v_points_[nn - 1]) / (alpha_k[nn] - alpha_k[nn - 1]) * std::exp(alpha_k[nn]);
            } else {
              g_k[nn] += (v_points_[nn] - v_points_[nn - 1]) / 2. * std::exp(alpha_k[nn]);
            }
          } // if (nn > 0)
          if (nn < dimRange - 1) {
            if (XT::Common::FloatCmp::ne(alpha_k[nn + 1], alpha_k[nn], taylor_tol_)) {
              g_k[nn] +=
                  (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                      * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                  - (v_points_[nn + 1] - v_points_[nn]) / (alpha_k[nn + 1] - alpha_k[nn]) * std::exp(alpha_k[nn]);
            } else {
              g_k[nn] += (v_points_[nn + 1] - v_points_[nn]) / 2. * std::exp(alpha_k[nn]);
            }
          } // if (nn < dimRange-1)
        } // nn
        g_k -= v;

        // calculate Hessian H
        H_k *= 0;
        for (size_t nn = 0; nn < dimRange; ++nn) {
          if (nn > 0) {
            if (XT::Common::FloatCmp::ne(alpha_k[nn], alpha_k[nn - 1], taylor_tol_)) {
              H_k[nn][nn - 1] =
                  (v_points_[nn] - v_points_[nn - 1])
                  * ((std::exp(alpha_k[nn]) + std::exp(alpha_k[nn - 1])) / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
                     - 2. * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                           / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3));
              H_k[nn][nn] +=
                  (v_points_[nn] - v_points_[nn - 1])
                  * ((-2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2) + 1. / (alpha_k[nn] - alpha_k[nn - 1]))
                         * std::exp(alpha_k[nn])
                     + 2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3)
                           * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1])));

            } else {
              H_k[nn][nn - 1] = (v_points_[nn] - v_points_[nn - 1]) / 6. * std::exp(alpha_k[nn]);
              H_k[nn][nn] += (v_points_[nn] - v_points_[nn - 1]) / 3. * std::exp(alpha_k[nn]);
            }
          } // if (nn > 0)
          if (nn < dimRange - 1) {
            if (XT::Common::FloatCmp::ne(alpha_k[nn + 1], alpha_k[nn], taylor_tol_)) {
              H_k[nn][nn + 1] =
                  (v_points_[nn + 1] - v_points_[nn])
                  * ((std::exp(alpha_k[nn + 1]) + std::exp(alpha_k[nn])) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                     - 2. * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                           / std::pow(alpha_k[nn + 1] - alpha_k[nn], 3));
              H_k[nn][nn] +=
                  (v_points_[nn + 1] - v_points_[nn])
                  * ((-2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2) - 1. / (alpha_k[nn + 1] - alpha_k[nn]))
                         * std::exp(alpha_k[nn])
                     + 2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 3)
                           * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn])));
            } else {
              H_k[nn][nn + 1] = (v_points_[nn + 1] - v_points_[nn]) / 6. * std::exp(alpha_k[nn]);
              H_k[nn][nn] += (v_points_[nn + 1] - v_points_[nn]) / 3. * std::exp(alpha_k[nn]);
            }
          } // if (nn < dimRange - 1)
        } // nn

        //        if (chol_flag == false && r == r_max)
        //          DUNE_THROW(Dune::NotImplemented, "Failure to converge!");
        // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
        //        if ((kk > k_0_ || chol_flag == false) && r < r_max)
        //          break;

        // calculate descent direction d_k;
        RangeType d_k(0), minus_g_k(g_k);
        minus_g_k *= -1;
        H_k.solve(d_k, minus_g_k);
        if (g_k.two_norm() < tau_ && std::exp(5 * d_k.one_norm()) < 1 + epsilon_gamma_) {
          alpha = alpha_k;
          goto outside_all_loops;
        } else {
          RangeFieldType zeta_k = 1;
          // backtracking line search
          while (zeta_k > epsilon_ * alpha_k.two_norm() / d_k.two_norm()) {

            // calculate alpha_new = alpha_k + zeta_k d_k
            auto alpha_new = d_k;
            alpha_new *= zeta_k;
            alpha_new += alpha_k;

            // calculate f(alpha_new)
            RangeFieldType f_new(0);
            for (size_t ii = 0; ii < dimRange - 1; ++ii) {
              if (XT::Common::FloatCmp::ne(alpha_new[ii + 1], alpha_new[ii], taylor_tol_))
                f_new += (v_points_[ii + 1] - v_points_[ii]) / (alpha_new[ii + 1] - alpha_new[ii])
                         * (std::exp(alpha_new[ii + 1]) - std::exp(alpha_new[ii]));
              else
                f_new += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_new[ii]);
            }
            f_new -= alpha_new * v;

            if (XT::Common::FloatCmp::le(f_new, f_k + xi_ * zeta_k * (g_k * d_k))) {
              alpha_k = alpha_new;
              f_k = f_new;
              break;
            }
            zeta_k = chi_ * zeta_k;
          } // backtracking linesearch while
        } // else (stopping conditions)
      } // k loop (Newton iterations)
    } // r loop (Regularization parameter)

    DUNE_THROW(NotImplemented, "Failed to converge");

  outside_all_loops:
    // store values as initial conditions for next time step on this entity
    alpha_cache_[id] = alpha;

    // TODO: calculate < \mu m G_\alpha(u) > !!!!!
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
          ret[nn] += -(v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                         * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
                     + (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
        } else {
          ret[nn] += (v_points_[nn] - v_points_[nn - 1]) / 2. * std::exp(alpha[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], taylor_tol_)) {
          ret[nn] += (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                         * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
                     - (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
        } else {
          ret[nn] += (v_points_[nn + 1] - v_points_[nn]) / 2. * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange-1)
    } // nn
    return ret;
  } // FluxRangeType evaluate(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  const typename GridType::LocalIdSet& id_set_;
  const RangeType v_points_;
  const RangeFieldType tau_;
  const RangeFieldType epsilon_gamma_;
  const RangeFieldType chi_;
  const RangeFieldType xi_;
  const std::vector<RangeFieldType> r_sequence_;
  const size_t k_0_;
  const size_t k_max_;
  const RangeFieldType epsilon_;
  const RangeFieldType taylor_tol_;
  const std::string name_;
  mutable std::unordered_map<typename GridType::LocalIdSet::IdType, RangeType> alpha_cache_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
