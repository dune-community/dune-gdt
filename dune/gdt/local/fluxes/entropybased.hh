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

namespace Dune {
namespace GDT {
namespace internal {


template <size_t dimRange>
Dune::FieldMatrix<double, dimRange, dimRange> unit_matrix()
{
  Dune::FieldMatrix<double, dimRange, dimRange> ret(0);
  for (size_t ii = 0; ii < dimRange; ++ii)
    ret[ii][ii] = 1.;
  return ret;
}

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
  for (int ii = int(A.N() - 1); ii >= 0; ii--) {
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


} // namespace internal


/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 */
template <class BasisfunctionType,
          class GridLayerType,
          class E,
          class D,
          size_t d,
          class U,
          class R,
          size_t rangeDim,
          size_t quadratureDim = d>
class EntropyBasedLocalFlux : public XT::Functions::LocalizableFluxFunctionInterface<E, D, d, U, 0, R, rangeDim, d>
{
  typedef XT::Functions::LocalizableFluxFunctionInterface<E, D, d, U, 0, R, rangeDim, d> BaseType;
  typedef EntropyBasedLocalFlux<BasisfunctionType, GridLayerType, E, D, d, U, R, rangeDim, quadratureDim> ThisType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::StateType;
  using typename BaseType::StateRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::PartialURangeType;
  using typename BaseType::LocalfunctionType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  static const size_t dimQuadrature = quadratureDim;
  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  using BasisValuesMatrixType = std::vector<StateRangeType>;
  typedef Dune::QuadratureRule<DomainFieldType, dimQuadrature> QuadratureRuleType;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& quadrature,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5},
      const size_t k_0 = 50,
      const size_t k_max = 200,
      const RangeFieldType epsilon = std::pow(2, -52),
      const MatrixType& T_minus_one = internal::unit_matrix<dimRange>(),
      const std::string name = static_id())
    : index_set_(grid_layer.indexSet())
    , basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , M_(quadrature_.size())
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
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      M_[ii] = basis_functions_.evaluate(quadrature_[ii].position());
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

    Localfunction(const EntityType& e,
                  const BasisfunctionType& basis_functions,
                  const QuadratureRuleType& quadrature,
                  const BasisValuesMatrixType& M,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const MatrixType& T_minus_one,
                  std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache,
                  std::unique_ptr<StateRangeType>& beta_cache,
                  std::unique_ptr<MatrixType>& T_cache,
                  std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_boundary,
                  std::unique_ptr<StateRangeType>& beta_cache_boundary,
                  std::unique_ptr<MatrixType>& T_cache_boundary)
      : LocalfunctionType(e)
      , basis_functions_(basis_functions)
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
      , alpha_cache_(alpha_cache)
      , beta_cache_(beta_cache)
      , T_cache_(T_cache)
      , alpha_cache_boundary_(alpha_cache_boundary)
      , beta_cache_boundary_(beta_cache_boundary)
      , T_cache_boundary_(T_cache_boundary)
    {
    }

    using LocalfunctionType::entity;

    StateRangeType
    get_alpha(const DomainType& /*x_local*/, const StateRangeType& u, const XT::Common::Parameter param) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      StateRangeType alpha;

      static size_t hitcounter = 0;
      // if value has already been calculated for these values, skip computation
      if (!boundary && alpha_cache_ && XT::Common::FloatCmp::eq(alpha_cache_->first, u)) {
        alpha = alpha_cache_->second;
        std::cout << hitcounter++ << std::endl;
      } else if (boundary && alpha_cache_boundary_ && XT::Common::FloatCmp::eq(alpha_cache_boundary_->first, u)) {
        alpha = alpha_cache_boundary_->second;
        std::cout << hitcounter++ << std::endl;
      } else {

        StateRangeType u_iso, alpha_iso;
        std::tie(u_iso, alpha_iso) = basis_functions_.calculate_isotropic_distribution(u);

        // define further variables
        bool chol_flag = false;
        StateRangeType g_k, beta_out;
        MatrixType T_k;
        BasisValuesMatrixType P_k(M_.size());

        const auto r_max = r_sequence_.back();
        for (const auto& r : r_sequence_) {
          StateRangeType beta_in = (beta_cache_ && !boundary)
                                       ? *beta_cache_
                                       : ((beta_cache_boundary_ && boundary) ? *beta_cache_boundary_ : alpha_iso);
          T_k = (T_cache_ && !boundary) ? *T_cache_
                                        : ((T_cache_boundary_ && boundary) ? *T_cache_boundary_ : T_minus_one_);
          // normalize u
          StateRangeType r_times_u_iso = u_iso;
          r_times_u_iso *= r;
          StateRangeType v = u;
          v *= 1 - r;
          v += r_times_u_iso;
          // calculate T_k u
          StateRangeType v_k;
          internal::solve_lower_triangular(T_k, v_k, v);
          // calculate values of basis p = T_k m
          for (size_t ii = 0; ii < M_.size(); ++ii)
            internal::solve_lower_triangular(T_k, P_k[ii], M_[ii]);
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
            StateRangeType error(0);
            for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
              auto m = M_[ll];
              StateRangeType Tinv_m(0);
              internal::solve_lower_triangular(T_k, Tinv_m, m);
              m *= std::exp(beta_out * Tinv_m) * quadrature_[ll].weight();
              error += m;
            }
            error -= v;
            // calculate descent direction d_k;
            StateRangeType d_k = g_k;
            d_k *= -1;
            StateRangeType T_k_inv_transp_d_k;
            try {
              internal::solve_lower_triangular_transposed(T_k, T_k_inv_transp_d_k, d_k);
            } catch (const Dune::FMatrixError& e) {
              if (r < r_max)
                break;
              else
                DUNE_THROW(Dune::FMatrixError, e.what());
            }
            if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.one_norm()) < 1 + epsilon_gamma_) {
              internal::solve_lower_triangular_transposed(T_k, alpha, beta_out);
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
        if (!boundary) {
          alpha_cache_ = std::make_unique<std::pair<StateRangeType, StateRangeType>>(std::make_pair(u, alpha));
          beta_cache_ = std::make_unique<StateRangeType>(beta_out);
          T_cache_ = std::make_unique<MatrixType>(T_k);
        } else {
          alpha_cache_boundary_ = std::make_unique<std::pair<StateRangeType, StateRangeType>>(std::make_pair(u, alpha));
          beta_cache_boundary_ = std::make_unique<StateRangeType>(beta_out);
          T_cache_boundary_ = std::make_unique<MatrixType>(T_k);
        }
      } // else ( value has not been calculated before )

      return alpha;
    }

    virtual size_t order(const XT::Common::Parameter& /*param*/) const override
    {
      return 1;
    }

    virtual void evaluate(const DomainType& x_local,
                          const StateRangeType& u,
                          RangeType& ret,
                          const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        const auto& position = quadrature_[ll].position();
        const auto& weight = quadrature_[ll].weight();
        const auto factor = std::exp(alpha * M_[ll]) * weight;
        for (size_t dd = 0; dd < dimDomain; ++dd)
          helper<dimDomain>::axpy(dd, ret, position[dd] * factor, M_[ll]);
      } // ll
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      // calculate ret[ii] = < omega[ii] m G_\alpha(u) >
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        const auto& position = quadrature_[ll].position();
        const auto& weight = quadrature_[ll].weight();
        const auto factor = std::exp(alpha * M_[ll]) * weight;
        ret.axpy(position[col] * factor, M_[ll]);
      } // ll
    } // void evaluate_col(...)


    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      MatrixType H_inv;
      calculate_hessian(alpha, M_, H_inv);
      H_inv.invert();
      helper<dimDomain>::partial_u(alpha, M_, H_inv, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      MatrixType H_inv;
      calculate_hessian(alpha, M_, H_inv);
      H_inv.invert();
      helper<dimDomain>::partial_u_col(col, alpha, M_, H_inv, ret, this);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    template <size_t domainDim = dimDomain, class anything = void>
    struct helper
    {
      static void axpy(const size_t dd, RangeType& ret, const RangeFieldType& factor, const StateRangeType& m)
      {
        for (size_t rr = 0; rr < dimRange; ++rr)
          ret[rr][dd] += factor * m[rr];
      }

      static void partial_u(const StateRangeType& alpha,
                            const BasisValuesMatrixType& M,
                            const MatrixType& H_inv,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        for (size_t dd = 0; dd < domainDim; ++dd) {
          entropy_flux->calculate_J(alpha, M, ret[dd], dd);
          ret[dd].rightmultiply(H_inv);
        }
      } // void partial_u(...)

      static void partial_u_col(const size_t col,
                                const StateRangeType& alpha,
                                const BasisValuesMatrixType& M,
                                const MatrixType& H_inv,
                                ColPartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(alpha, M, ret, col);
        ret.rightmultiply(H_inv);
      } // void partial_u_col(...)
    }; // class helper<...>

    template <class anything>
    struct helper<1, anything>
    {
      static void axpy(const size_t dd, RangeType& ret, const StateRangeType& m, const RangeFieldType& factor)
      {
        assert(dd == 0);
        ret.axpy(factor, m);
      }

      static void partial_u(const StateRangeType& alpha,
                            const BasisValuesMatrixType& M,
                            const MatrixType& H_inv,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(alpha, M, ret, 0);
        ret.rightmultiply(H_inv);
      } // void partial_u(...)

      static void partial_u_col(const size_t col,
                                const StateRangeType& alpha,
                                const BasisValuesMatrixType& M,
                                const MatrixType& H_inv,
                                PartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        assert(col == 0);
        partial_u(alpha, M, H_inv, ret, entropy_flux);
      } // void partial_u(...)
    }; // class helper<1, ...>

    void calculate_hessian(const StateRangeType& alpha, const BasisValuesMatrixType& M, MatrixType& H) const
    {
      std::fill(H.begin(), H.end(), 0);
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        auto m_times_factor = M[ll];
        m_times_factor *= std::exp(alpha * M[ll]) * quadrature_[ll].weight();
        // add m m^T * factor
        for (size_t rr = 0; rr < dimRange; ++rr)
          for (size_t cc = rr; cc < dimRange; ++cc)
            H[rr][cc] += M[ll][rr] * m_times_factor[cc];
      } // quadrature points for loop
      // symmetric update for lower triangular part of H
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < rr; ++cc)
          H[rr][cc] = H[cc][rr];
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    void calculate_J(const StateRangeType& alpha,
                     const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<double, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimRangeCols);
      std::fill(J_dd.begin(), J_dd.end(), 0);
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        const auto& v = quadrature_[ll].position();
        auto m_times_factor = M[ll];
        m_times_factor *= v[dd] * std::exp(alpha * M[ll]) * quadrature_[ll].weight();
        // add m m^T * factor
        for (size_t rr = 0; rr < dimRange; ++rr)
          for (size_t cc = rr; cc < dimRange; ++cc)
            J_dd[rr][cc] += M[ll][rr] * m_times_factor[cc];
      } // quadrature points for loop
      // symmetric update for lower triangular part of J
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < rr; ++cc)
          J_dd[rr][cc] = J_dd[cc][rr];
    } // void calculate_J(...)

    void change_basis(bool& chol_flag,
                      const StateRangeType& beta_in,
                      StateRangeType& v_k,
                      BasisValuesMatrixType& P_k,
                      MatrixType& T_k,
                      StateRangeType& g_k,
                      StateRangeType& beta_out) const
    {
      MatrixType H(0), L(0);
      calculate_hessian(beta_in, P_k, H);
      chol_flag = cholesky_L(H, L);
      if (chol_flag == false)
        return;
      const auto P_k_copy = P_k;
      const auto v_k_copy = v_k;
      for (size_t ll = 0; ll < P_k.size(); ++ll)
        internal::solve_lower_triangular(L, P_k[ll], P_k_copy[ll]);
      T_k.rightmultiply(L);
      L.mtv(beta_in, beta_out);
      internal::solve_lower_triangular(L, v_k, v_k_copy);
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

    const BasisfunctionType& basis_functions_;
    const QuadratureRuleType& quadrature_;
    const BasisValuesMatrixType& M_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const MatrixType& T_minus_one_;
    const std::string name_;
    // constructor)
    std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_;
    std::unique_ptr<StateRangeType>& beta_cache_;
    std::unique_ptr<MatrixType>& T_cache_;
    std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_boundary_;
    std::unique_ptr<StateRangeType>& beta_cache_boundary_;
    std::unique_ptr<MatrixType>& T_cache_boundary_;
  }; // class Localfunction>

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           basis_functions_,
                                           quadrature_,
                                           M_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           T_minus_one_,
                                           alpha_cache_[index],
                                           beta_cache_[index],
                                           T_cache_[index],
                                           alpha_cache_[index_set_.size(0) + index],
                                           beta_cache_[index_set_.size(0) + index],
                                           T_cache_[index_set_.size(0) + index]);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  StateRangeType evaluate_kinetic_integral(const EntityType& entity,
                                           const DomainType& x_local_entity,
                                           const StateRangeType& u_i,
                                           const EntityType& neighbor,
                                           const DomainType& x_local_neighbor,
                                           const StateRangeType u_j,
                                           const DomainType& n_ij,
                                           const XT::Common::Parameter& param,
                                           const XT::Common::Parameter& param_neighbor) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = local_function(entity);
    const auto local_function_neighbor = local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param);
    const auto alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor);
    StateRangeType ret(0);
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& position = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      const auto& m = M_[ll];
      const auto factor = position * n_ij > 0 ? std::exp(alpha_i * m) * weight : std::exp(alpha_j * m) * weight;
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        auto contribution = m;
        contribution *= position[dd] * factor * n_ij[dd];
        ret += contribution;
      }
    }
    return ret;
  } // StateRangeType evaluate_kinetic_integral(...)

  const typename GridLayerType::IndexSet& index_set_;
  const BasisfunctionType& basis_functions_;
  const QuadratureRuleType quadrature_;
  BasisValuesMatrixType M_;
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
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as much matrices or vectors as needed (see
  // constructor)
  mutable std::vector<std::unique_ptr<std::pair<StateRangeType, StateRangeType>>> alpha_cache_;
  mutable std::vector<std::unique_ptr<StateRangeType>> beta_cache_;
  mutable std::vector<std::unique_ptr<MatrixType>> T_cache_;
};

#if 0
/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 * domainDim, rangeDim, rangeDimCols are the respective dimensions of pde solution u, not the dimensions of \mathbf{f}.
 */
template <class GridLayerType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class EntropyBasedLocalFluxHatFunctions1D : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedLocalFluxHatFunctions1D<GridLayerType, E, D, d, R, rangeDim, rC> ThisType;

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
      const GridLayerType& grid_layer,
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
    : global_index_set_(grid_layer, 0)
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
  const Dune::GlobalIndexSet<GridLayerType> global_index_set_;
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
template <class GridLayerType, class E, class D, size_t d, class R, size_t rangeDim, size_t rC>
class EntropyBasedLocalFluxHatFunctions3d : public AnalyticalFluxInterface<E, D, d, R, rangeDim, rC>
{
  typedef AnalyticalFluxInterface<E, D, d, R, rangeDim, rC> BaseType;
  typedef EntropyBasedLocalFluxHatFunctions3d<GridLayerType, E, D, d, R, rangeDim, rC> ThisType;

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
      const GridLayerType& grid_layer,
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
    : global_index_set_(grid_layer, 0)
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

  const Dune::GlobalIndexSet<GridLayerType> global_index_set_;
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
