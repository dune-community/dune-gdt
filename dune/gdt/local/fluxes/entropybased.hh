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

#include <dune/grid/utility/globalindexset.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/math.hh>

#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/test/hyperbolic/problems/fokkerplanck/pointsource.hh>

namespace Dune {
namespace GDT {


template <size_t dimRange>
Dune::FieldMatrix<double, dimRange, dimRange> unit_matrix()
{
  Dune::FieldMatrix<double, dimRange, dimRange> ret(0);
  for (size_t ii = 0; ii < dimRange; ++ii)
    ret[ii][ii] = 1.;
  return ret;
}

template <class MatrixType>
MatrixType sparse_unit_matrix(const size_t rows, const size_t cols)
{
  XT::LA::SparsityPatternDefault pattern(rows);
  for (size_t ii = 0; ii < rows; ++ii)
    pattern.insert(ii, ii);
  MatrixType ret(rows, cols, pattern);
  for (size_t ii = 0; ii < rows; ++ii)
    XT::Common::MatrixAbstraction<MatrixType>::set_entry(ret, ii, ii, 1.);
  return ret;
}

template <class EntityType,
          class DomainType,
          class RangeType,
          class FluxRangeType,
          class QuadratureRuleType,
          class BasisValuesMatrixType,
          size_t dimDomain>
struct EntropyBasedLocalFluxEvaluator
{
  static FluxRangeType evaluate(const RangeType& u,
                                const EntityType& /*entity*/,
                                const DomainType& /*x_local*/,
                                const double /*t*/,
                                const QuadratureRuleType& quadrature,
                                const RangeType& alpha,
                                const BasisValuesMatrixType& M)
  {
    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    FluxRangeType ret(0);
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto& position = quadrature[ll].position();
      const auto& mu = position[0];
      const auto& phi = position[1];
      const auto& weight = quadrature[ll].weight();
      std::vector<double> omega(3);
      omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
      omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
      omega[2] = mu;
      RangeType m(0);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        m = M[ll];
        m *= omega[dd] * std::exp(alpha * m) * weight;
        for (size_t rr = 0; rr < u.size(); ++rr)
          ret[rr][dd] += m[rr];
      }
    }
    return ret;
  } // FluxRangeType evaluate(...)
};

template <class EntityType,
          class DomainType,
          class RangeType,
          class FluxRangeType,
          class QuadratureRuleType,
          class BasisValuesMatrixType>
struct EntropyBasedLocalFluxEvaluator<EntityType,
                                      DomainType,
                                      RangeType,
                                      FluxRangeType,
                                      QuadratureRuleType,
                                      BasisValuesMatrixType,
                                      1>
{
  static FluxRangeType evaluate(const RangeType& /*u*/,
                                const EntityType& /*entity*/,
                                const DomainType& /*x_local*/,
                                const double /*t*/,
                                const QuadratureRuleType& quadrature,
                                const RangeType& alpha,
                                const BasisValuesMatrixType& M)
  {


    // calculate < \mu m G_\alpha(u) >
    FluxRangeType ret(0);
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto& position = quadrature[ll].position();
      const auto& mu = position[0];
      const auto& weight = quadrature[ll].weight();
      auto m = M[ll];
      m *= mu * std::exp(alpha * m) * weight;
      ret += m;
    }
    return ret;
  }
};

template <size_t dimRange, size_t dimDomain>
struct FluxRangeTypeConverter
{
  static Dune::FieldMatrix<double, dimRange, dimDomain>
  convert(const Dune::FieldMatrix<double, dimRange, dimDomain>& ret)
  {
    return ret;
  }
};

template <size_t dimRange>
struct FluxRangeTypeConverter<dimRange, 1>
{
  static Dune::FieldVector<double, dimRange> convert(const Dune::FieldMatrix<double, dimRange, 1>& in)
  {
    Dune::FieldVector<double, dimRange> ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = in[ii][0];
    return ret;
  }
};

template <size_t dimRange, size_t dimDomain>
struct FluxJacobianRangeTypeConverter
{
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> FluxJacobianRangeType;
  static FluxJacobianRangeType convert(const FluxJacobianRangeType& ret)
  {
    return ret;
  }
};

template <size_t dimRange>
struct FluxJacobianRangeTypeConverter<dimRange, 1>
{
  typedef typename Dune::FieldMatrix<double, dimRange, dimRange> FluxJacobianRangeType;
  static FluxJacobianRangeType convert(const Dune::FieldVector<FluxJacobianRangeType, 1>& in)
  {
    return in[0];
  }
};

template <class MatrixType, class VectorType>
void solve_lower_triangular(const MatrixType& A, VectorType& x, const VectorType& b)
{
  VectorType& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // forward solve
  for (size_t ii = 0; ii < A.N(); ++ii) {
    for (size_t jj = 0; jj < ii; ++jj)
      rhs[ii] -= A[ii][jj] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
}

template <class MatrixType, class VectorType>
void solve_lower_triangular_transposed(const MatrixType& A, VectorType& x, const VectorType& b)
{
  VectorType& rhs = x; // use x to store rhs
  rhs = b; // copy data
  // backsolve
  double min_eigval(std::abs(A[0][0]));
  double max_eigval = min_eigval;
  for (int ii = A.N() - 1; ii >= 0; ii--) {
    auto abs = std::abs(A[ii][ii]);
    min_eigval = std::min(abs, min_eigval);
    max_eigval = std::max(abs, max_eigval);
    for (size_t jj = ii + 1; jj < A.N(); jj++)
      rhs[ii] -= A[jj][ii] * x[jj];
    x[ii] = rhs[ii] / A[ii][ii];
  }
  if (XT::Common::FloatCmp::eq(min_eigval, 0.) || max_eigval / min_eigval > 1e10)
    DUNE_THROW(Dune::FMatrixError, "A is singular!");
}

#if 1
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of
 \mathbf{f}.
 */
template <class GridViewType,
          class E,
          class D,
          size_t d,
          class R,
          size_t rangeDim,
          size_t rC,
          bool quadrature_is_cartesian = true>
class EntropyBasedLocalFlux : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedLocalFlux<GridViewType, E, D, d, R, rangeDim, rC, quadrature_is_cartesian> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using BasisValuesMatrixType = std::vector<FieldVector<RangeFieldType, dimRange>>;
  using typename BaseType::FluxRangeType;
  using typename BaseType::FluxJacobianRangeType;
  typedef Dune::QuadratureRule<DomainFieldType,
                               quadrature_is_cartesian ? dimDomain : (dimDomain == 1 ? 1 : dimDomain - 1)>
      QuadratureRuleType;
  typedef std::function<std::pair<RangeType, RangeType>(const RangeType&)> IsotropicDistributionCalculatorType;

  explicit EntropyBasedLocalFlux(
      const GridViewType& grid_view,
      const QuadratureRuleType& quadrature,
      const BasisValuesMatrixType& M,
      const IsotropicDistributionCalculatorType isotropic_dist_calculator,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = unit_matrix<dimRange>(),
      const std::string name = static_id())
    : index_set_(grid_view.indexSet())
    , quadrature_(quadrature)
    , M_(M)
    , isotropic_dist_calculator_(isotropic_dist_calculator)
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
    , alpha_cache_(2 * index_set_.size(0))
    , beta_cache_(2 * index_set_.size(0))
    , T_cache_(2 * index_set_.size(0))
  {
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);
    MatrixType H;
    Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> ret;
    calculate_hessian(alpha, M_, H);
    calculate_J(alpha, M_, ret);
    H.invert();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      ret[dd].rightmultiply(H);
    return ret;
  }

  RangeType get_alpha(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    // get initial multiplier and basis matrix from last time step
    // in the numerical flux, we are setting x_local to DomainType(200) if we are actually not on the entity,
    // but on the boundary (YaspGrid does not support ghost entities, thus we use this hack)
    const auto index = index_set_.index(entity) + index_set_.size(0) * (x_local[0] > 100);
    RangeType alpha;

    // if value has already been calculated for this entity at this time, skip computation
    //    if (alpha_cache_[index] && XT::Common::FloatCmp::eq(alpha_cache_[index]->first, t)) {
    //      alpha = alpha_cache_[index]->second;
    //    } else {

    RangeType u_iso, alpha_iso;
    std::tie(u_iso, alpha_iso) = isotropic_dist_calculator_(u);

    // define further variables
    bool chol_flag = false;
    RangeType g_k, beta_out;
    MatrixType T_k;
    BasisValuesMatrixType P_k(M_.size());

    const auto r_max = r_sequence_.back();
    for (const auto& r : r_sequence_) {
      RangeType beta_in = beta_cache_[index] ? *(beta_cache_[index]) : alpha_iso;
      T_k = T_cache_[index] ? *(T_cache_[index]) : T_minus_one_;
      // normalize u
      RangeType r_times_u_iso = u_iso;
      r_times_u_iso *= r;
      RangeType v = u;
      v *= 1 - r;
      v += r_times_u_iso;
      // calculate T_k u
      RangeType v_k;
      solve_lower_triangular(T_k, v_k, v);
      // calculate values of basis p = T_k m
      for (size_t ii = 0; ii < M_.size(); ++ii)
        solve_lower_triangular(T_k, P_k[ii], M_[ii]);
      // calculate f_0
      RangeFieldType f_k(0);
      for (size_t ll = 0; ll < quadrature_.size(); ++ll)
        f_k += quadrature_[ll].weight() * std::exp(beta_in * P_k[ll]);
      f_k -= beta_in * v_k;

      for (size_t kk = 0; kk < k_max_; ++kk) {
        change_basis(chol_flag, beta_in, v_k, P_k, T_k, g_k, beta_out);
        if (chol_flag == false && r == r_max)
          DUNE_THROW(Dune::MathError, "Failure to converge!");
        // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
        if ((kk > k_0_ || chol_flag == false) && r < r_max)
          break;
        // calculate current error
        RangeType error(0);
        for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
          auto m = M_[ll];
          RangeType Tinv_m(0);
          solve_lower_triangular(T_k, Tinv_m, m);
          m *= std::exp(beta_out * Tinv_m) * quadrature_[ll].weight();
          error += m;
        }
        error -= v;
        // calculate descent direction d_k;
        RangeType d_k = g_k;
        d_k *= -1;
        RangeType T_k_inv_transp_d_k;
        try {
          solve_lower_triangular_transposed(T_k, T_k_inv_transp_d_k, d_k);
        } catch (const Dune::FMatrixError& e) {
          if (r < r_max)
            break;
          else
            DUNE_THROW(Dune::FMatrixError, e.what());
        }
        if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.one_norm()) < 1 + epsilon_gamma_) {
          solve_lower_triangular_transposed(T_k, alpha, beta_out);
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

    DUNE_THROW(MathError, "Failed to converge");

  outside_all_loops:
    // store values as initial conditions for next time step on this entity
    alpha_cache_[index] = std::make_unique<std::pair<double, RangeType>>(std::make_pair(t, alpha));
    beta_cache_[index] = std::make_unique<RangeType>(beta_out);
    T_cache_[index] = std::make_unique<MatrixType>(T_k);
    //    } // else ( value has not been calculated before )

    return alpha;
  }

  virtual FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);
    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> ret(0);
    Dune::FieldVector<RangeFieldType, dimDomain> omega;
    RangeType m;
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& position = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      if (quadrature_is_cartesian) {
        for (size_t dd = 0; dd < dimDomain; ++dd)
          omega[dd] = position[dd];
      } else {
        const auto& mu = position[0];
        const auto& phi = position[1];
        omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
        omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
        omega[2] = mu;
      }
      const auto factor = std::exp(alpha * m) * weight;
      m = M_[ll];
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        m *= omega[dd] * factor;
        for (size_t rr = 0; rr < dimRange; ++rr)
          ret[rr][dd] += m[rr];
      }
    }
    return FluxRangeTypeConverter<dimRange, dimDomain>::convert(ret);
    //    return ret;
  } // FluxRangeType evaluate(...)

  virtual RangeType calculate_flux_integral(const RangeType& u_i,
                                            const E& entity,
                                            const DomainType& x_local_entity,
                                            const RangeType u_j,
                                            const E& neighbor,
                                            const DomainType& x_local_neighbor,
                                            const DomainType& n_ij,
                                            const double t) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, entity, x_local_entity, t);
    const auto alpha_j = get_alpha(u_j, neighbor, x_local_neighbor, t);
    RangeType ret(0);
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& position = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      Dune::FieldVector<double, dimDomain> omega;
      if (quadrature_is_cartesian) {
        for (size_t dd = 0; dd < dimDomain; ++dd)
          omega[dd] = position[dd];
      } else {
        const auto& mu = position[0];
        const auto& phi = position[1];
        omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
        omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
        omega[2] = mu;
      }
      const auto& m = M_[ll];
      const auto factor = position * n_ij > 0 ? std::exp(alpha_i * m) * weight : std::exp(alpha_j * m) * weight;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        auto contribution = m;
        contribution *= omega[dd] * factor * n_ij[dd];
        ret += contribution;
      }
    }
    return ret;
  } // RangeType calculate_flux_integral(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  void calculate_hessian(const RangeType& alpha, const BasisValuesMatrixType& M, MatrixType& H) const
  {
    H *= 0;
    MatrixType tmp;
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& m = M[ll];
      auto m_times_factor = M[ll];
      m_times_factor *= std::exp(alpha * m) * quadrature_[ll].weight();
      // calculate p p^T
      for (size_t ii = 0; ii < dimRange; ++ii) {
        tmp[ii] = m_times_factor;
        tmp[ii] *= m[ii];
      }
      H += tmp;
    } // quadrature points for loop
  }

  // J = df/dalpha is the derivative of the flux with respect to alpha.
  // As F = (f_1, f_2, f_3) is matrix-valued
  // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
  // vector-valued),
  // the derivative is vector of matrices (df_1/dalpha, df_2/dalpha, ...)
  void calculate_J(const RangeType& alpha,
                   const BasisValuesMatrixType& M,
                   Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain>& J) const
  {
    std::fill(J.begin(), J.end(), MatrixType(0));
    MatrixType tmp;
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& m = M[ll];
      const auto& v = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        auto m_times_factor = M[ll];
        m_times_factor *= v[dd] * std::exp(alpha * m) * weight;
        // calculate p p^T
        for (size_t ii = 0; ii < dimRange; ++ii) {
          tmp[ii] = m_times_factor;
          tmp[ii] *= m[ii];
        }
        J[dd] += tmp;
      } // dd
    } // quadrature points for loop
  }

  void change_basis(bool& chol_flag,
                    const RangeType& beta_in,
                    RangeType& v_k,
                    BasisValuesMatrixType& P_k,
                    MatrixType& T_k,
                    RangeType& g_k,
                    RangeType& beta_out) const
  {
    MatrixType H(0), L(0);
    calculate_hessian(beta_in, P_k, H);
    chol_flag = cholesky_L(H, L);
    if (chol_flag == false)
      return;
    const auto P_k_copy = P_k;
    const auto v_k_copy = v_k;
    for (size_t ll = 0; ll < P_k.size(); ++ll)
      solve_lower_triangular(L, P_k[ll], P_k_copy[ll]);
    T_k.rightmultiply(L);
    L.mtv(beta_in, beta_out);
    solve_lower_triangular(L, v_k, v_k_copy);
    g_k = v_k;
    g_k *= -1;
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      // assumes that the first basis function is constant
      //      g_k[0] += std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight() * P_k[0][0];
      auto contribution = P_k[ll];
      contribution *= std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight();
      g_k += contribution;
    }
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

      if (XT::Common::FloatCmp::le(xDiag, RangeFieldType(0)))
        return false;

      rii = std::sqrt(xDiag);

      RangeFieldType invrii = RangeFieldType(1) / rii;
      for (int ll = ii + 1; ll < size; ++ll) {
        RangeFieldType x = H[ll][ii];
        for (int jj = 0; jj < ii; ++jj)
          x -= L[ii][jj] * L[ll][jj];
        L[ll][ii] = invrii * x;
      }
    }
    return true;
  }

  const typename GridViewType::IndexSet& index_set_;
  const QuadratureRuleType quadrature_;
  const BasisValuesMatrixType M_;
  const IsotropicDistributionCalculatorType isotropic_dist_calculator_;
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
  // TODO: concurrent writes to different locations of std::vector are thread-safe in most implementations of
  // std::vector,
  // but that is not guaranteed by the standard.
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as much matrices or vectors as needed (see
  // constructor)
  mutable std::vector<std::unique_ptr<std::pair<double, RangeType>>> alpha_cache_;
  mutable std::vector<std::unique_ptr<RangeType>> beta_cache_;
  mutable std::vector<std::unique_ptr<MatrixType>> T_cache_;
};
#endif

#if 0
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of
 \mathbf{f}.
 */
template <class GridViewType,
          class E,
          class D,
          size_t d,
          class R,
          size_t rangeDim,
          size_t rC,
          bool quadrature_is_cartesian = true>
class AdaptiveEntropyBasedLocalFlux : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef AdaptiveEntropyBasedLocalFlux<GridViewType, E, D, d, R, rangeDim, rC, quadrature_is_cartesian> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using BasisValuesMatrixType = std::vector<FieldVector<RangeFieldType, dimRange>>;
  using typename BaseType::FluxRangeType;
  using typename BaseType::FluxJacobianRangeType;
  typedef Dune::QuadratureRule<DomainFieldType,
                               quadrature_is_cartesian ? dimDomain : (dimDomain == 1 ? 1 : dimDomain - 1)>
      QuadratureRuleType;
  typedef std::function<std::pair<RangeType, RangeType>(const RangeType&)> IsotropicDistributionCalculatorType;

  explicit AdaptiveEntropyBasedLocalFlux(
      const GridViewType& grid_view,
      const std::vector<QuadratureRuleType>& quadratures,
      const BasisValuesMatrixType& M,
      const IsotropicDistributionCalculatorType isotropic_dist_calculator,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = unit_matrix<dimRange>(),
      const std::string name = static_id())
    : index_set_(grid_view.indexSet())
    , quadratures_(quadratures)
    , M_(M)
    , isotropic_dist_calculator_(isotropic_dist_calculator)
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
    , alpha_cache_(2 * index_set_.size(0))
    , beta_cache_(2 * index_set_.size(0))
    , T_cache_(2 * index_set_.size(0))
  {
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& /*u*/, const E& /*entity*/, const DomainType& /*x_local*/, const double /*t*/ = 0) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  RangeType get_alpha(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    // get initial multiplier and basis matrix from last time step
    // in the numerical flux, we are setting x_local to DomainType(200) if we are actually not on the entity,
    // but on the boundary (YaspGrid does not support ghost entities, thus we use this hack)
    const auto index = index_set_.index(entity) + index_set_.size(0) * (x_local[0] > 100);
    RangeType alpha;

    // if value has already been calculated for this entity at this time, skip computation
    if (alpha_cache_[index] && XT::Common::FloatCmp::eq(alpha_cache_[index]->first, t)) {
      alpha = alpha_cache_[index]->second;
    } else {

      RangeType u_iso, alpha_iso;
      std::tie(u_iso, alpha_iso) = isotropic_dist_calculator_(u);

      // define further variables
      bool chol_flag = false;
      RangeType g_k, beta_out;
      MatrixType T_k;
      size_t current_quadrature_level(2);
      size_t num_quadratures(quadratures_.size());
      size_t next_quadrature_level =
          current_quadrature_level == num_quadratures - 1 ? current_quadrature_level : current_quadrature_level + 1;
      const auto* current_quadrature = &(quadratures_[current_quadrature_level]);
      const auto* next_quadrature = &(quadratures_[next_quadrature_level]);
      BasisValuesMatrixType P_k(next_quadrature->size());

      const auto r_max = r_sequence_.back();
      for (const auto& r : r_sequence_) {
        RangeType beta_in = beta_cache_[index] ? *(beta_cache_[index]) : alpha_iso;
        T_k = T_cache_[index] ? *(T_cache_[index]) : T_minus_one_;
        // normalize u
        RangeType r_times_u_iso = u_iso;
        r_times_u_iso *= r;
        RangeType v = u;
        v *= 1 - r;
        v += r_times_u_iso;
        // calculate T_k u
        RangeType v_k;
        solve_lower_triangular(T_k, v_k, v);
        // calculate values of basis p = T_k m
        // also for finer quadrature, as we will need it in the end anyway
        for (size_t ii = 0; ii < next_quadrature->size(); ++ii)
          solve_lower_triangular(T_k, P_k[ii], M_[ii]);
        // calculate f_0
        RangeFieldType f_k(0);
        for (size_t ll = 0; ll < current_quadrature->size(); ++ll)
          f_k += (*current_quadrature)[ll].weight() * std::exp(beta_in * P_k[ll]);
        f_k -= beta_in * v_k;

        std::vector<RangeType> beta_iterates;
        for (size_t kk = 0; kk < k_max_; ++kk) {
          beta_iterates.emplace_back(beta_in);
          change_basis(chol_flag, beta_in, v_k, P_k, T_k, g_k, beta_out, current_quadrature);
          if (chol_flag == false && r == r_max)
            DUNE_THROW(Dune::MathError, "Failure to converge!");
          // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
          if ((kk > k_0_ || chol_flag == false) && r < r_max)
            break;
          // calculate current error
          RangeType error(0);
          for (size_t ll = 0; ll < current_quadrature->size(); ++ll) {
            auto m = M_[ll];
            RangeType Tinv_m(0);
            solve_lower_triangular(T_k, Tinv_m, m);
            m *= std::exp(beta_out * Tinv_m) * (*current_quadrature)[ll].weight();
            error += m;
          }
          error -= v;
          // calculate descent direction d_k;
          RangeType d_k = g_k;
          d_k *= -1;
          RangeType T_k_inv_transp_d_k;
          try {
            solve_lower_triangular_transposed(T_k, T_k_inv_transp_d_k, d_k);
          } catch (const Dune::FMatrixError& e) {
            if (r < r_max)
              break;
            else
              DUNE_THROW(Dune::FMatrixError, e.what());
          }
          if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.one_norm()) < 1 + epsilon_gamma_) {
            solve_lower_triangular_transposed(T_k, alpha, beta_out);
            goto outside_all_loops;
          } else {
            RangeFieldType zeta_k = 1;
            beta_in = beta_out;
            // backtracking line search
            bool same_quadrature = true;
            while (zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm()) {
              RangeFieldType f(0);
              auto beta_new = d_k;
              beta_new *= zeta_k;
              beta_new += beta_out;
              for (size_t ll = 0; ll < current_quadrature->size(); ++ll)
                f += (*current_quadrature)[ll].weight() * std::exp(beta_new * P_k[ll]);
              f -= beta_new * v_k;
              if (XT::Common::FloatCmp::le(f, f_k + xi_ * zeta_k * (g_k * d_k))) {

                // check with finer quadrature if max quadrature level is not reached
                if (next_quadrature->size() > current_quadrature->size()) {
                  RangeFieldType f_new_quad(0), f_k_new_quad(0);
                  //                  RangeType g_k_new_quad = v_k;
                  //                  g_k_new_quad *= -1;
                  for (size_t ll = 0; ll < next_quadrature->size(); ++ll) {
                    f_new_quad += (*next_quadrature)[ll].weight() * std::exp(beta_new * P_k[ll]);
                    f_k_new_quad += (*next_quadrature)[ll].weight() * std::exp(beta_out * P_k[ll]);
                    //                    auto g_k_contribution = P_k[ll];
                    //                    g_k_contribution *= std::exp(beta_out * P_k[ll]) *
                    //                    (*current_quadrature)[ll].weight();
                    //                    g_k_new_quad += g_k_contribution;
                  }
                  f_new_quad -= beta_new * v_k;
                  f_k_new_quad -= beta_out * v_k;

                  if (XT::Common::FloatCmp::le(f_new_quad, f_k_new_quad + xi_ * zeta_k * (g_k * d_k))) {

                  } else {
                    // continue linesearch with finer quadrature
                    f_k = f_k_new_quad;
                    refine_quadrature(
                        current_quadrature_level, next_quadrature_level, current_quadrature, next_quadrature, T_k, P_k);
                    zeta_k = chi_ * zeta_k;
                    same_quadrature = false;
                    continue;
                  }
                }

                beta_in = beta_new;
                f_k = f;
                break;
              }
              zeta_k = chi_ * zeta_k;
              if (XT::Common::FloatCmp::le(zeta_k, epsilon_ * beta_out.two_norm() / d_k.two_norm())) {
                if (same_quadrature) {
                  if (beta_iterates.size() == 0) {
                    break;
                  }
                  beta_in = beta_iterates.back();
                  beta_iterates.pop_back();
                }
                change_basis(chol_flag, beta_in, v_k, P_k, T_k, g_k, beta_out, current_quadrature);
                if (chol_flag == false && beta_iterates.size() == 0 && r == r_max)
                  DUNE_THROW(Dune::MathError, "Failed to converge nnnn!");
                // exit inner for loop to increase r if to many iterations are used or cholesky decomposition fails
                if (chol_flag == false) {
                  beta_in = beta_iterates.back();
                  beta_iterates.pop_back();
                  change_basis(chol_flag, beta_in, v_k, P_k, T_k, g_k, beta_out, current_quadrature);
                  d_k = g_k;
                  d_k *= -1.;
                  zeta_k = 1.;
                  continue;
                }
                d_k = g_k;
                d_k *= -1.;
                zeta_k = 1.;
                same_quadrature = true;
              }
            } // backtracking linesearch while
          } // else (stopping conditions)
        } // k loop (Newton iterations)
      } // r loop (Regularization parameter)

      DUNE_THROW(MathError, "Failed to converge");

    outside_all_loops:
      // store values as initial conditions for next time step on this entity
      alpha_cache_[index] = std::make_unique<std::pair<double, RangeType>>(std::make_pair(t, alpha));
      beta_cache_[index] = std::make_unique<RangeType>(beta_out);
      T_cache_[index] = std::make_unique<MatrixType>(T_k);
    } // else ( value has not been calculated before )

    return alpha;
  }

  void refine_quadrature(size_t& current_quadrature_level,
                         size_t& next_quadrature_level,
                         const QuadratureRuleType* current_quadrature,
                         const QuadratureRuleType* next_quadrature,
                         const MatrixType& T_k,
                         BasisValuesMatrixType& P_k) const
  {
    if (next_quadrature_level == quadratures_.size() - 1)
      next_quadrature_level = current_quadrature_level;
    current_quadrature = &(quadratures_[++current_quadrature_level]);
    next_quadrature = &(quadratures_[++next_quadrature_level]);
    P_k.resize(next_quadrature->size());
    if (current_quadrature->size() != next_quadrature->size())
      for (size_t ii = current_quadrature->size(); ii < next_quadrature->size(); ++ii)
        solve_lower_triangular(T_k, P_k[ii], M_[ii]);
  }

  virtual FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);
    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> ret(0);
    Dune::FieldVector<RangeFieldType, dimDomain> omega;
    RangeType m;
    for (size_t ll = 0; ll < quadratures_.back().size(); ++ll) {
      const auto& position = quadratures_.back()[ll].position();
      const auto& weight = quadratures_.back()[ll].weight();
      if (quadrature_is_cartesian) {
        for (size_t dd = 0; dd < dimDomain; ++dd)
          omega[dd] = position[dd];
      } else {
        const auto& mu = position[0];
        const auto& phi = position[1];
        omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
        omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
        omega[2] = mu;
      }
      const auto factor = std::exp(alpha * m) * weight;
      m = M_[ll];
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        m *= omega[dd] * factor;
        for (size_t rr = 0; rr < dimRange; ++rr)
          ret[rr][dd] += m[rr];
      }
    }
    return FluxRangeTypeConverter<dimRange, dimDomain>::convert(ret);
    //    return ret;
  } // FluxRangeType evaluate(...)

  virtual RangeType calculate_flux_integral(const RangeType& u_i,
                                            const E& entity,
                                            const DomainType& x_local_entity,
                                            const RangeType u_j,
                                            const E& neighbor,
                                            const DomainType& x_local_neighbor,
                                            const DomainType& n_ij,
                                            const double t) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, entity, x_local_entity, t);
    const auto alpha_j = get_alpha(u_j, neighbor, x_local_neighbor, t);
    RangeType ret(0);
    const auto& quadrature = quadratures_.back();
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto& position = quadrature[ll].position();
      const auto& weight = quadrature[ll].weight();
      Dune::FieldVector<double, dimDomain> omega;
      if (quadrature_is_cartesian) {
        for (size_t dd = 0; dd < dimDomain; ++dd)
          omega[dd] = position[dd];
      } else {
        const auto& mu = position[0];
        const auto& phi = position[1];
        omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
        omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
        omega[2] = mu;
      }
      const auto& m = M_[ll];
      const auto factor = position * n_ij > 0 ? std::exp(alpha_i * m) * weight : std::exp(alpha_j * m) * weight;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        auto contribution = m;
        contribution *= omega[dd] * factor * n_ij[dd];
        ret += contribution;
      }
    }
    return ret;
  } // RangeType calculate_flux_integral(...)

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
                    RangeType& beta_out,
                    const QuadratureRuleType* current_quadrature) const
  {
    MatrixType H(0), tmp(0), L(0);
    for (size_t ll = 0; ll < current_quadrature->size(); ++ll) {
      const auto& p = P_k[ll];
      auto p_times_factor = P_k[ll];
      p_times_factor *= std::exp(beta_in * p) * (*current_quadrature)[ll].weight();
      // calculate p p^T
      for (size_t ii = 0; ii < dimRange; ++ii) {
        tmp[ii] = p_times_factor;
        tmp[ii] *= p[ii];
      }
      H += tmp;
    } // quadrature points for loop
    chol_flag = cholesky_L(H, L);
    if (chol_flag == false)
      return;
    const auto P_k_copy = P_k;
    const auto v_k_copy = v_k;
    for (size_t ll = 0; ll < P_k.size(); ++ll)
      solve_lower_triangular(L, P_k[ll], P_k_copy[ll]);
    T_k.rightmultiply(L);
    L.mtv(beta_in, beta_out);
    solve_lower_triangular(L, v_k, v_k_copy);
    g_k = v_k;
    g_k *= -1;
    for (size_t ll = 0; ll < current_quadrature->size(); ++ll) {
      // assumes that the first basis function is constant
      //      g_k[0] += std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight() * P_k[0][0];
      auto contribution = P_k[ll];
      contribution *= std::exp(beta_out * P_k[ll]) * (*current_quadrature)[ll].weight();
      g_k += contribution;
    }
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

      if (XT::Common::FloatCmp::le(xDiag, RangeFieldType(0)))
        return false;

      rii = std::sqrt(xDiag);

      RangeFieldType invrii = RangeFieldType(1) / rii;
      for (int ll = ii + 1; ll < size; ++ll) {
        RangeFieldType x = H[ll][ii];
        for (int jj = 0; jj < ii; ++jj)
          x -= L[ii][jj] * L[ll][jj];
        L[ll][ii] = invrii * x;
      }
    }
    return true;
  }

  const typename GridViewType::IndexSet& index_set_;
  const std::vector<QuadratureRuleType>& quadratures_;
  const BasisValuesMatrixType& M_;
  const IsotropicDistributionCalculatorType isotropic_dist_calculator_;
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
  // TODO: concurrent writes to different locations of std::vector are thread-safe in most implementations of
  // std::vector,
  // but that is not guaranteed by the standard.
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as much matrices or vectors as needed (see
  // constructor)
  mutable std::vector<std::unique_ptr<std::pair<double, RangeType>>> alpha_cache_;
  mutable std::vector<std::unique_ptr<RangeType>> beta_cache_;
  mutable std::vector<std::unique_ptr<MatrixType>> T_cache_;
};
#endif


#if 0
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of
 \mathbf{f}.
 */
template <class GridViewType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class AdaptiveEntropyBasedLocalFlux : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef AdaptiveEntropyBasedLocalFlux<GridViewType, E, D, d, R, rangeDim, rC> ThisType;

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
  typedef std::function<std::pair<RangeType, RangeType>(const RangeType&)> IsotropicDistributionCalculatorType;
//  typedef typename GDT::Hyperbolic::Problems::AdaptiveQuadrature<DomainType, RangeType, RangeType>
//      AdaptiveQuadratureType;
//  typedef typename AdaptiveQuadratureType::QuadraturePointType QuadraturePointType;
  typedef typename GridViewType::IndexSet IndexSetType;

  explicit AdaptiveEntropyBasedLocalFlux(
      const GridViewType& grid_view,
      const IsotropicDistributionCalculatorType isotropic_dist_calculator,
      AdaptiveQuadratureType& adaptive_quadrature,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = unit_matrix<dimRange>(),
      const std::string name = static_id())
    : index_set_(grid_view.indexSet())
    , isotropic_dist_calculator_(isotropic_dist_calculator)
    , adaptive_quadrature_(adaptive_quadrature)
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
    , alpha_cache_(2 * index_set_.size(0))
    , beta_cache_(2 * index_set_.size(0))
    , T_cache_(2 * index_set_.size(0))
  {
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& /*u*/, const E& /*entity*/, const DomainType& /*x_local*/, const double /*t*/ = 0) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  RangeType get_alpha(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    // get initial multiplier and basis matrix from last time step
    // in the numerical flux, we are setting x_local to DomainType(200) if we are actually not on the entity,
    // but on the boundary (YaspGrid does not support ghost entities, thus we use this hack)
    const auto index = index_set_.index(entity) + index_set_.size(0) * (x_local[0] > 100);
    RangeType alpha;

    // if value has already been calculated for this entity at this time, skip computation
    if (alpha_cache_[index] && XT::Common::FloatCmp::eq(alpha_cache_[index]->first, t)) {
      alpha = alpha_cache_[index]->second;
    } else {

      RangeType u_iso, alpha_iso;
      std::tie(u_iso, alpha_iso) = isotropic_dist_calculator_(u);

      // define further variables
      bool chol_flag = false;
      RangeType g_k, beta_out;
      MatrixType T_k;

      const auto r_max = r_sequence_.back();
      for (const auto& r : r_sequence_) {
        adaptive_quadrature_.reset();
        RangeType beta_in = beta_cache_[index] ? *(beta_cache_[index]) : alpha_iso;
        T_k = T_cache_[index] ? *(T_cache_[index]) : T_minus_one_;
        // normalize u
        RangeType r_times_u_iso = u_iso;
        r_times_u_iso *= r;
        RangeType v = u;
        v *= 1 - r;
        v += r_times_u_iso;
        // calculate T_k u
        RangeType v_k;
        T_k.solve(v_k, v);

        // calculate f_0
        RangeFieldType f_k;
        try {
          f_k = adaptive_quadrature_.template calculate_integral<RangeFieldType>(
              [&](const QuadraturePointType& quadpoint, const RangeType& /*m*/, const RangeType& p) {
                return std::exp(beta_in * p) * quadpoint.weight();
              },
              [&](const QuadraturePointType& /*quadpoint*/, const RangeType& m) {
                RangeType ret(0);
                solve_lower_triangular(T_k, ret, m);
                return ret;
              },
              true);
        } catch (const Dune::MathError& e) {
          if (r < r_max)
            continue;
          else
            DUNE_THROW(Dune::MathError, e.what());
        }
        f_k -= beta_in * v_k;

        for (size_t kk = 0; kk < k_max_; ++kk) {
          change_basis(chol_flag, beta_in, v_k, T_k, g_k, beta_out);
          if (chol_flag == false && r == r_max)
            DUNE_THROW(Dune::MathError, "Failure to converge!");
          // exit inner for loop to increase r if too many iterations are used or cholesky decomposition fails
          if ((kk > k_0_ || chol_flag == false) && r < r_max)
            break;
          // calculate current error
          RangeType error;
          try {
            error = adaptive_quadrature_.template calculate_integral<RangeType>(
                [&](const QuadraturePointType& quadpoint, const RangeType& m, const RangeType& /*p*/) {
                  RangeType Tinv_m(0);
                  solve_lower_triangular(T_k, Tinv_m, m);
                  auto ret = m;
                  ret *= std::exp(beta_out * Tinv_m) * quadpoint.weight();
                  return ret;
                },
                [&](const QuadraturePointType& /*quadpoint*/, const RangeType& m) {
                  ;
                  RangeType ret(0);
                  solve_lower_triangular(T_k, ret, m);
                  return ret;
                },
                false);
          } catch (const Dune::MathError& e) {
            if (r < r_max)
              break;
            else
              DUNE_THROW(Dune::MathError, "");
          }
          error -= v;
          // calculate descent direction d_k;
          RangeType d_k = g_k;
          d_k *= -1;
          RangeType T_k_inv_transp_d_k;
          try {
            solve_lower_triangular_transposed(T_k, T_k_inv_transp_d_k, d_k);
          } catch (const Dune::FMatrixError& e) {
            if (r < r_max)
              break;
            else
              DUNE_THROW(Dune::FMatrixError, e.what());
          }
          if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.one_norm()) < 1 + epsilon_gamma_) {
            solve_lower_triangular_transposed(T_k, alpha, beta_out);
            goto outside_all_loops;
          } else {
            RangeFieldType zeta_k = 1;
            beta_in = beta_out;
            // backtracking line search
            while (zeta_k > epsilon_ * beta_out.two_norm() / d_k.two_norm()) {
              auto beta_new = d_k;
              beta_new *= zeta_k;
              beta_new += beta_out;

              RangeFieldType f;
              try {
                f = adaptive_quadrature_.template calculate_integral<RangeFieldType>(
                    [&](const QuadraturePointType& quadpoint, const RangeType& /*m*/, const RangeType& p) {
                      return std::exp(beta_new * p) * quadpoint.weight();
                    },
                    [&](const QuadraturePointType& /*quadpoint*/, const RangeType& m) {
                      RangeType ret(0);
                      solve_lower_triangular(T_k, ret, m);
                      return ret;
                    },
                    false);
              } catch (const Dune::MathError& e) {
                if (r < r_max) {
                  kk = 10000;
                  break;
                } else
                  DUNE_THROW(Dune::MathError, "");
              }
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

      DUNE_THROW(MathError, "Failed to converge");

    outside_all_loops:
      // store values as initial conditions for next time step on this entity
      alpha_cache_[index] = std::make_unique<std::pair<double, RangeType>>(std::make_pair(t, alpha));
      beta_cache_[index] = std::make_unique<RangeType>(beta_out);
      T_cache_[index] = std::make_unique<MatrixType>(T_k);
    } // else ( value has not been calculated before )

    return alpha;
  }

  virtual FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    //    const auto alpha = get_alpha(u, entity, x_local, t);
    //    // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
    //    Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> ret(0);
    //    Dune::FieldVector<RangeFieldType, dimDomain> omega;
    //    RangeType m;
    //    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
    //      const auto& position = quadrature_[ll].position();
    //      const auto& weight = quadrature_[ll].weight();
    //      if (quadrature_is_cartesian) {
    //        for (size_t dd = 0; dd < dimDomain; ++dd)
    //          omega[dd] = position[dd];
    //      } else {
    //        const auto& mu = position[0];
    //        const auto& phi = position[1];
    //        omega[0] = std::sqrt(1. - mu * mu) * std::cos(phi);
    //        omega[1] = std::sqrt(1. - mu * mu) * std::sin(phi);
    //        omega[2] = mu;
    //      }
    //      const auto factor = std::exp(alpha * m) * weight;
    //      m = M_[ll];
    //      for (size_t dd = 0; dd < dimDomain; ++dd) {
    //        m *= omega[dd] * factor;
    //        for (size_t rr = 0; rr < dimRange; ++rr)
    //          ret[rr][dd] += m[rr];
    //      }
    //    }

    // return FluxRangeTypeConverter<dimRange, dimDomain>::convert(ret);
    return FluxRangeType(0);

  } // FluxRangeType evaluate(...)

  virtual RangeType calculate_flux_integral(const RangeType& u_i,
                                            const E& entity,
                                            const DomainType& x_local_entity,
                                            const RangeType u_j,
                                            const E& neighbor,
                                            const DomainType& x_local_neighbor,
                                            const DomainType& n_ij,
                                            const double t) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto alpha_i = get_alpha(u_i, entity, x_local_entity, t);
    const auto alpha_j = get_alpha(u_j, neighbor, x_local_neighbor, t);
    return adaptive_quadrature_.template calculate_integral<RangeType>(
        [&](const QuadraturePointType& quadpoint, const RangeType& m, const RangeType& /*p*/) {
          const auto factor = quadpoint.position() * n_ij > 0 ? std::exp(alpha_i * m) * quadpoint.weight()
                                                              : std::exp(alpha_j * m) * quadpoint.weight();
          RangeType ret(0);
          for (size_t dd = 0; dd < dimDomain; ++dd) {
            auto summand = m;
            summand *= quadpoint.position()[dd] * factor * n_ij[dd];
            ret += summand;
          } // dd
          return ret;
        },
        [&](const QuadraturePointType& /*quadpoint*/, const RangeType& /*m*/) { return RangeType(0); },
        false);
  } // RangeType calculate_flux_integral(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  void change_basis(bool& chol_flag,
                    const RangeType& beta_in,
                    RangeType& v_k,
                    MatrixType& T_k,
                    RangeType& g_k,
                    RangeType& beta_out) const
  {
    MatrixType tmp, L(0), H;

    try {
      H = adaptive_quadrature_.template calculate_integral<MatrixType>(
          [&](const QuadraturePointType& quadpoint, const RangeType& /*m*/, const RangeType& p) {
            auto p_scaled = p;
            p_scaled *= std::exp(beta_in * p) * quadpoint.weight();
            for (size_t ii = 0; ii < dimRange; ++ii) {
              tmp[ii] = p_scaled;
              tmp[ii] *= p[ii];
            }
            return tmp;
          },
          [&](const QuadraturePointType& /*quadpoint*/, const RangeType& m) {
            RangeType ret(0);
            solve_lower_triangular(T_k, ret, m);
            return ret;
          },
          false);
    } catch (const Dune::MathError&) {
      chol_flag = false;
      return;
    }

    chol_flag = cholesky_L(H, L);
    if (chol_flag == false)
      return;
    const auto v_k_copy = v_k;
    T_k.rightmultiply(L);
    L.mtv(beta_in, beta_out);
    solve_lower_triangular(L, v_k, v_k_copy);

    try {
      g_k = adaptive_quadrature_.template calculate_integral<RangeType>(
          [&](const QuadraturePointType& quadpoint, const RangeType& /*m*/, const RangeType& p) {
            auto p_scaled = p;
            p_scaled *= std::exp(beta_out * p) * quadpoint.weight();
            return p_scaled;
          },
          [&](const QuadraturePointType& /*quadpoint*/, const RangeType& m) {
            RangeType ret(0);
            solve_lower_triangular(T_k, ret, m);
            return ret;
          },
          true);
      g_k -= v_k;
    } catch (const Dune::MathError&) {
      chol_flag = false;
      return;
    }
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

      if (XT::Common::FloatCmp::le(xDiag, RangeFieldType(0)))
        return false;

      rii = std::sqrt(xDiag);

      RangeFieldType invrii = RangeFieldType(1) / rii;
      for (int ll = ii + 1; ll < size; ++ll) {
        RangeFieldType x = H[ll][ii];
        for (int jj = 0; jj < ii; ++jj)
          x -= L[ii][jj] * L[ll][jj];
        L[ll][ii] = invrii * x;
      }
    }
    return true;
  }

  const IndexSetType& index_set_;
  const IsotropicDistributionCalculatorType isotropic_dist_calculator_;
  AdaptiveQuadratureType& adaptive_quadrature_;
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
  mutable std::vector<std::unique_ptr<std::pair<double, RangeType>>> alpha_cache_;
  mutable std::vector<std::unique_ptr<RangeType>> beta_cache_;
  mutable std::vector<std::unique_ptr<MatrixType>> T_cache_;
};
#endif

#if 0
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of \mathbf{f}.
 */
template <class GridViewType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class EntropyBasedLocalFluxHatFunctions1D : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedLocalFluxHatFunctions1D<GridViewType, E, D, d, R, rangeDim, rC> ThisType;

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

  explicit EntropyBasedLocalFluxHatFunctions1D(
      const GridViewType& grid_view,
      const RangeType v_points,
      const RangeFieldType tau = 1e-7,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const RangeFieldType taylor_tol = 1e-4,
      const size_t taylor_order = 10,
      const std::string name = static_id())
    : global_index_set_(grid_view, 0)
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
    , taylor_order_(taylor_order)
    , name_(name)
    , alpha_cache_(2 * global_index_set_.size(0))
  {
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& /*u*/, const E& /*entity*/, const DomainType& /*x_local*/, const double /*t*/ = 1) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  RangeType get_alpha(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    // in the numerical flux, we are setting x_local to DomainType(200) if we are actually not on the entity,
    // but on the boundary (YaspGrid does not support ghost entities, thus we use this hack)
    const auto index = global_index_set_.index(entity) + global_index_set_.size(0) * (x_local[0] > 100);
    RangeType alpha;

    // if value has already been calculated for this entity at this time, skip computation
    if (alpha_cache_[index] && alpha_cache_[index]->first == t) {
      alpha = alpha_cache_[index]->second;
    } else {
      // get initial multiplier and basis matrix from last time step
      RangeType alpha_iso(1);
      RangeFieldType psi_iso(0);
      for (size_t ii = 0; ii < dimRange; ++ii)
        psi_iso += u[ii];
      psi_iso /= 2.;
      alpha_iso *= std::log(psi_iso);

      // define further variables
      RangeType g_k;
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
        RangeType alpha_k = alpha_cache_[index] ? alpha_cache_[index]->second : alpha_iso;
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
          else {
            RangeFieldType taylorsum = 0.;
            for (size_t ll = 1; ll <= taylor_order_; ++ll)
              taylorsum += std::pow(alpha_k[ii + 1] - alpha_k[ii], ll - 1.) / XT::Common::factorial(ll);
            f_k += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_k[ii]) * taylorsum;
          }
        }
        f_k -= alpha_k * v;

        for (size_t kk = 0; kk < k_max_; ++kk) {
          // exit inner for loop to increase r if too many iterations are used
          if (kk > k_0_ && r < r_max)
            break;

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
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 2.) / XT::Common::factorial(ll);
                g_k[nn] += taylorsum * std::exp(alpha_k[nn]) * (v_points_[nn] - v_points_[nn - 1]);
              }
            } // if (nn > 0)
            if (nn < dimRange - 1) {
              if (XT::Common::FloatCmp::ne(alpha_k[nn + 1], alpha_k[nn], taylor_tol_)) {
                g_k[nn] +=
                    (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                        * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                    - (v_points_[nn + 1] - v_points_[nn]) / (alpha_k[nn + 1] - alpha_k[nn]) * std::exp(alpha_k[nn]);
              } else {
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 2.) / XT::Common::factorial(ll);
                g_k[nn] += taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
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
                H_k[nn][nn] =
                    (v_points_[nn] - v_points_[nn - 1])
                    * ((-2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2) + 1. / (alpha_k[nn] - alpha_k[nn - 1]))
                           * std::exp(alpha_k[nn])
                       + 2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3)
                             * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1])));

              } else {
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 2.)
                               * (1. / XT::Common::factorial(ll) - 2. / XT::Common::factorial(ll + 1));
                H_k[nn][nn - 1] = taylorsum * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
                taylorsum = 0.;
                for (size_t ll = 3; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 3.) * 2. / XT::Common::factorial(ll);
                H_k[nn][nn] = taylorsum * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
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
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 2.)
                               * (1. / XT::Common::factorial(ll) - 2. / XT::Common::factorial(ll + 1));
                H_k[nn][nn + 1] = taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
                taylorsum = 0.;
                for (size_t ll = 3; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 3.) * 2. / XT::Common::factorial(ll);
                H_k[nn][nn] += taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
              }
            } // if (nn < dimRange - 1)
          } // nn

          // calculate descent direction d_k;
          RangeType d_k(0), minus_g_k(g_k);
          minus_g_k *= -1;
          try {
            H_k.solve(d_k, minus_g_k);
          } catch (const Dune::FMatrixError& err) {
            if (r < r_max) {
              break;
            } else {
              DUNE_THROW(Dune::FMatrixError, "Failure to converge!");
            }
          }

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
                else {
                  RangeFieldType taylorsum = 0.;
                  for (size_t ll = 1; ll <= taylor_order_; ++ll)
                    taylorsum += std::pow(alpha_new[ii + 1] - alpha_new[ii], ll - 1.) / XT::Common::factorial(ll);
                  f_new += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_new[ii]) * taylorsum;
                }
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

      DUNE_THROW(MathError, "Failed to converge");

    outside_all_loops:
      // store values as initial conditions for next time step on this entity
      alpha_cache_[index] = std::make_unique<std::pair<double, RangeType>>(std::make_pair(t, alpha));
    }
    return alpha;
  }

  // integrals can be evaluated exactly for hatfunctions
  FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);

    // calculate < \mu m G_\alpha(u) >
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
          ret[nn] +=
              2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                  * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
              + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                    * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn] - v_points_[nn - 1], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn - 1] * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], taylor_tol_)) {
          ret[nn] +=
              -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                  * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
              + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                    * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn + 1] - v_points_[nn], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn + 1] * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn

    return ret;
  } // FluxRangeType evaluate(...)

  virtual FluxRangeType calculate_flux_integral(const RangeType& u_i,
                                                const E& entity,
                                                const DomainType& x_local_entity,
                                                const RangeType u_j,
                                                const E& neighbor,
                                                const DomainType& x_local_neighbor,
                                                const DomainType& n_ij,
                                                const double t) const
  {
    assert(v_points_.size() % 2 && "Not implemented for odd number of points!");
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto alpha_i = get_alpha(u_i, entity, x_local_entity, t);
    const auto alpha_j = get_alpha(u_j, neighbor, x_local_neighbor, t);
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        const auto& alpha =
            XT::Common::FloatCmp::ge(n_ij[0] * (v_points_[nn - 1] + v_points_[nn]) / 2., 0.) ? alpha_i : alpha_j;
        if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
          ret[nn] +=
              2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                  * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
              + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                    * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn] - v_points_[nn - 1], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn - 1] * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        const auto& alpha =
            XT::Common::FloatCmp::ge(n_ij[0] * (v_points_[nn] + v_points_[nn + 1]) / 2., 0.) ? alpha_i : alpha_j;
        if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], taylor_tol_)) {
          ret[nn] +=
              -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                  * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
              + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                    * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn + 1] * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn

    ret *= n_ij[0];

    return ret;
  } // FluxRangeType calculate_flux_integral(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  const Dune::GlobalIndexSet<GridViewType> global_index_set_;
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
  const size_t taylor_order_;
  const std::string name_;
  mutable std::vector<std::unique_ptr<std::pair<double, RangeType>>> alpha_cache_;
};


/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of \mathbf{f}.
 */
template <class GridViewType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class EntropyBasedLocalFluxHatFunctions3d : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedLocalFluxHatFunctions3d<GridViewType, E, D, d, R, rangeDim, rC> ThisType;

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
  typedef Dune::GDT::Hyperbolic::Problems::SphericalTriangulation<RangeFieldType> TriangulationType;

  explicit EntropyBasedLocalFluxHatFunctions3d(
      const GridViewType& grid_view,
      const TriangulationType triangulation,
      const RangeFieldType tau = 1e-7,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const size_t taylor_order = 10,
      const std::string name = static_id())
    : global_index_set_(grid_view, 0)
    , triangulation_(triangulation)
    , tau_(tau)
    , epsilon_gamma_(epsilon_gamma)
    , chi_(chi)
    , xi_(xi)
    , r_sequence_(r_sequence)
    , k_0_(k_0)
    , k_max_(k_max)
    , epsilon_(epsilon)
    , taylor_order_(taylor_order)
    , name_(name)
    , alpha_cache_(2 * global_index_set_.size(0))
  {
  }

  virtual FluxJacobianRangeType
  jacobian(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);
  }

  RangeType get_alpha(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    // in the numerical flux, we are setting x_local to DomainType(200) if we are actually not on the entity,
    // but on the boundary (YaspGrid does not support ghost entities, thus we use this hack)
    const auto index = global_index_set_.index(entity) + global_index_set_.size(0) * (x_local[0] > 100);
    RangeType alpha;

    // if value has already been calculated for this entity at this time, skip computation
    if (alpha_cache_[index] && alpha_cache_[index]->first == t) {
      alpha = alpha_cache_[index]->second;
    } else {
      // get initial multiplier and basis matrix from last time step
      RangeType alpha_iso(1);
      RangeFieldType psi_iso(0);
      for (size_t ii = 0; ii < dimRange; ++ii)
        psi_iso += u[ii];
      psi_iso /= 4. * M_PI;
      alpha_iso *= std::log(psi_iso);

      // define further variables
      RangeType g_k;
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
        RangeType alpha_k = alpha_cache_[index] ? alpha_cache_[index]->second : alpha_iso;
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
          else {
            RangeFieldType taylorsum = 0.;
            for (size_t ll = 1; ll <= taylor_order_; ++ll)
              taylorsum += std::pow(alpha_k[ii + 1] - alpha_k[ii], ll - 1.) / XT::Common::factorial(ll);
            f_k += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_k[ii]) * taylorsum;
          }
        }
        f_k -= alpha_k * v;

        for (size_t kk = 0; kk < k_max_; ++kk) {
          // exit inner for loop to increase r if too many iterations are used
          if (kk > k_0_ && r < r_max)
            break;

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
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 2.) / XT::Common::factorial(ll);
                g_k[nn] += taylorsum * std::exp(alpha_k[nn]) * (v_points_[nn] - v_points_[nn - 1]);
              }
            } // if (nn > 0)
            if (nn < dimRange - 1) {
              if (XT::Common::FloatCmp::ne(alpha_k[nn + 1], alpha_k[nn], taylor_tol_)) {
                g_k[nn] +=
                    (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                        * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                    - (v_points_[nn + 1] - v_points_[nn]) / (alpha_k[nn + 1] - alpha_k[nn]) * std::exp(alpha_k[nn]);
              } else {
                RangeFieldType taylorsum = 0.;
                for (size_t ll = 2; ll <= taylor_order_; ++ll)
                  taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 2.) / XT::Common::factorial(ll);
                g_k[nn] += taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
              }
            } // if (nn < dimRange-1)
          } // nn
          g_k -= v;

          // calculate Hessian H
          calculate_hessian(alpha_k, H_k);

          // calculate descent direction d_k;
          RangeType d_k(0), minus_g_k(g_k);
          minus_g_k *= -1;
          try {
            H_k.solve(d_k, minus_g_k);
          } catch (const Dune::FMatrixError& err) {
            if (r < r_max) {
              break;
            } else {
              DUNE_THROW(Dune::FMatrixError, "Failure to converge!");
            }
          }

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
                else {
                  RangeFieldType taylorsum = 0.;
                  for (size_t ll = 1; ll <= taylor_order_; ++ll)
                    taylorsum += std::pow(alpha_new[ii + 1] - alpha_new[ii], ll - 1.) / XT::Common::factorial(ll);
                  f_new += (v_points_[ii + 1] - v_points_[ii]) * std::exp(alpha_new[ii]) * taylorsum;
                }
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

      DUNE_THROW(MathError, "Failed to converge");

    outside_all_loops:
      // store values as initial conditions for next time step on this entity
      alpha_cache_[index] = std::make_unique<std::pair<double, RangeType>>(std::make_pair(t, alpha));
    }
    return alpha;
  }

  // integrals can be evaluated exactly for hatfunctions
  FluxRangeType evaluate(const RangeType& u, const E& entity, const DomainType& x_local, const double t) const
  {
    const auto alpha = get_alpha(u, entity, x_local, t);

    // calculate < \mu m G_\alpha(u) >
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
          ret[nn] +=
              2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                  * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
              + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                    * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn] - v_points_[nn - 1], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn - 1] * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], taylor_tol_)) {
          ret[nn] +=
              -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                  * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
              + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                    * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn + 1] - v_points_[nn], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn + 1] * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn

    return ret;
  } // FluxRangeType evaluate(...)

  virtual FluxRangeType calculate_flux_integral(const RangeType& u_i,
                                                const E& entity,
                                                const DomainType& x_local_entity,
                                                const RangeType u_j,
                                                const E& neighbor,
                                                const DomainType& x_local_neighbor,
                                                const DomainType& n_ij,
                                                const double t) const
  {
    assert(v_points_.size() % 2 && "Not implemented for odd number of points!");
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto alpha_i = get_alpha(u_i, entity, x_local_entity, t);
    const auto alpha_j = get_alpha(u_j, neighbor, x_local_neighbor, t);
    RangeType ret(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        const auto& alpha =
            XT::Common::FloatCmp::ge(n_ij[0] * (v_points_[nn - 1] + v_points_[nn]) / 2., 0.) ? alpha_i : alpha_j;
        if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
          ret[nn] +=
              2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                  * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
              + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                    * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * 2 * std::pow(v_points_[nn] - v_points_[nn - 1], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn - 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn - 1] * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        const auto& alpha =
            XT::Common::FloatCmp::ge(n_ij[0] * (v_points_[nn] + v_points_[nn + 1]) / 2., 0.) ? alpha_i : alpha_j;
        if (XT::Common::FloatCmp::ne(alpha[nn + 1], alpha[nn], taylor_tol_)) {
          ret[nn] +=
              -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                  * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
              + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                    * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                       - 2 * v_points_[nn] * std::exp(alpha[nn]))
              - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn]) * std::exp(alpha[nn]);
        } else {
          RangeFieldType taylorsum = 0;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 3.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) * std::exp(alpha[nn]);
          taylorsum = 0;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha[nn + 1] - alpha[nn], ll - 2.) / XT::Common::factorial(ll);
          ret[nn] += taylorsum * v_points_[nn + 1] * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn

    ret *= n_ij[0];

    return ret;
  } // FluxRangeType calculate_flux_integral(...)

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

private:
  void calculate_hessian(const RangeType& alpha_k, MatrixType& H_k)
  {
    H_k *= 0;
    for (size_t nn = 0; nn < dimRange; ++nn) {
      if (nn > 0) {
        if (XT::Common::FloatCmp::ne(alpha_k[nn], alpha_k[nn - 1], taylor_tol_)) {
          H_k[nn][nn - 1] =
              (v_points_[nn] - v_points_[nn - 1])
              * ((std::exp(alpha_k[nn]) + std::exp(alpha_k[nn - 1])) / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
                 - 2. * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1]))
                       / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3));
          H_k[nn][nn] = (v_points_[nn] - v_points_[nn - 1])
                        * ((-2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2) + 1. / (alpha_k[nn] - alpha_k[nn - 1]))
                               * std::exp(alpha_k[nn])
                           + 2. / std::pow(alpha_k[nn] - alpha_k[nn - 1], 3)
                                 * (std::exp(alpha_k[nn]) - std::exp(alpha_k[nn - 1])));

        } else {
          RangeFieldType taylorsum = 0.;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 2.)
                         * (1. / XT::Common::factorial(ll) - 2. / XT::Common::factorial(ll + 1));
          H_k[nn][nn - 1] = taylorsum * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
          taylorsum = 0.;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha_k[nn - 1] - alpha_k[nn], ll - 3.) * 2. / XT::Common::factorial(ll);
          H_k[nn][nn] = taylorsum * (v_points_[nn] - v_points_[nn - 1]) * std::exp(alpha_k[nn]);
        }
      } // if (nn > 0)
      if (nn < dimRange - 1) {
        if (XT::Common::FloatCmp::ne(alpha_k[nn + 1], alpha_k[nn], taylor_tol_)) {
          H_k[nn][nn + 1] =
              (v_points_[nn + 1] - v_points_[nn])
              * ((std::exp(alpha_k[nn + 1]) + std::exp(alpha_k[nn])) / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
                 - 2. * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn]))
                       / std::pow(alpha_k[nn + 1] - alpha_k[nn], 3));
          H_k[nn][nn] += (v_points_[nn + 1] - v_points_[nn])
                         * ((-2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2) - 1. / (alpha_k[nn + 1] - alpha_k[nn]))
                                * std::exp(alpha_k[nn])
                            + 2. / std::pow(alpha_k[nn + 1] - alpha_k[nn], 3)
                                  * (std::exp(alpha_k[nn + 1]) - std::exp(alpha_k[nn])));
        } else {
          RangeFieldType taylorsum = 0.;
          for (size_t ll = 2; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 2.)
                         * (1. / XT::Common::factorial(ll) - 2. / XT::Common::factorial(ll + 1));
          H_k[nn][nn + 1] = taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
          taylorsum = 0.;
          for (size_t ll = 3; ll <= taylor_order_; ++ll)
            taylorsum += std::pow(alpha_k[nn + 1] - alpha_k[nn], ll - 3.) * 2. / XT::Common::factorial(ll);
          H_k[nn][nn] += taylorsum * (v_points_[nn + 1] - v_points_[nn]) * std::exp(alpha_k[nn]);
        }
      } // if (nn < dimRange - 1)
    } // nn
  } // void calculate_hessian(...)

  const Dune::GlobalIndexSet<GridViewType> global_index_set_;
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
  const size_t taylor_order_;
  const std::string name_;
  mutable std::vector<std::unique_ptr<std::pair<double, RangeType>>> alpha_cache_;
};
#endif


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
