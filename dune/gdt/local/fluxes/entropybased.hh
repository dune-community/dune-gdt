// This file is part of the dune-gdt project:
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

#include <dune/xt/la/container/algorithms.hh>
#include <dune/xt/la/container/common.hh>

#include <dune/gdt/local/fluxes/interfaces.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/hatfunctions.hh>

namespace Dune {
namespace GDT {


template <class MatrixType, size_t size>
static const std::unique_ptr<const MatrixType> unit_matrix = XT::LA::get_unit_matrix<MatrixType>(size, 0);

/** Analytical flux \mathbf{f}(\mathbf{u}) = < \mu \mathbf{m} G_{\hat{\alpha}(\mathbf{u})} >,
 * for the notation see
 * Alldredge, Hauck, O'Leary, Tits, "Adaptive change of basis in entropy-based moment closures for linear kinetic
 * equations"
 */
template <class BasisfunctionImp, class GridLayerImp, class U, size_t quadratureDim = BasisfunctionImp::dimDomain>
class EntropyBasedLocalFlux
    : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                             typename BasisfunctionImp::DomainFieldType,
                                                             BasisfunctionImp::dimFlux,
                                                             U,
                                                             0,
                                                             typename BasisfunctionImp::RangeFieldType,
                                                             BasisfunctionImp::dimRange,
                                                             BasisfunctionImp::dimFlux>
{
  typedef typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                                   typename BasisfunctionImp::DomainFieldType,
                                                                   BasisfunctionImp::dimFlux,
                                                                   U,
                                                                   0,
                                                                   typename BasisfunctionImp::RangeFieldType,
                                                                   BasisfunctionImp::dimRange,
                                                                   BasisfunctionImp::dimFlux>
      BaseType;
  typedef EntropyBasedLocalFlux ThisType;

public:
  typedef BasisfunctionImp BasisfunctionType;
  typedef GridLayerImp GridLayerType;
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
  using FieldMatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  //  using MatrixType = FieldMatrix<RangeFieldType, dimRange, dimRange>;
  //  using MatrixType = XT::LA::CommonSparseMatrix<RangeFieldType>;
  using MatrixType = XT::LA::CommonSparseMatrixCsc<RangeFieldType>;
  //  using SparseMatrixType = XT::LA::CommonSparseMatrix<RangeFieldType>;
  typedef typename XT::LA::CommonSparseVector<RangeFieldType> VectorType;
  using BasisValuesMatrixType = std::vector<VectorType>;
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
      const MatrixType& T_minus_one = *unit_matrix<MatrixType, dimRange>,
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
    , mutexes_(index_set_.size(0))
  {
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      M_[ii] = VectorType(basis_functions_.evaluate(quadrature_[ii].position()), true, size_t(0));
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
                  std::unique_ptr<std::pair<StateRangeType, VectorType>>& alpha_cache,
                  std::unique_ptr<VectorType>& beta_cache,
                  std::unique_ptr<MatrixType>& T_cache,
                  std::unique_ptr<std::pair<StateRangeType, VectorType>>& alpha_cache_boundary,
                  std::unique_ptr<VectorType>& beta_cache_boundary,
                  std::unique_ptr<MatrixType>& T_cache_boundary,
                  std::mutex& mutex)
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
      , mutex_(mutex)
    {
    }

    using LocalfunctionType::entity;

    VectorType
    get_alpha(const DomainType& /*x_local*/, const StateRangeType& u, const XT::Common::Parameter param) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      // get initial multiplier and basis matrix from last time step
      static thread_local VectorType alpha(dimRange, size_t(0));

      static std::atomic<size_t> hitcounter(0);
      // if value has already been calculated for these values, skip computation
      mutex_.lock();
      if (!boundary && alpha_cache_ && XT::Common::FloatCmp::eq(alpha_cache_->first, u)) {
        alpha.deep_copy(alpha_cache_->second);
        //        std::cout << hitcounter++ << std::endl;
        mutex_.unlock();
        return alpha;
      } else if (boundary && alpha_cache_boundary_ && XT::Common::FloatCmp::eq(alpha_cache_boundary_->first, u)) {
        alpha.deep_copy(alpha_cache_boundary_->second);
        //        std::cout << hitcounter++ << std::endl;
        mutex_.unlock();
        return alpha;
      } else {

        StateRangeType u_iso, alpha_iso;
        std::tie(u_iso, alpha_iso) = basis_functions_.calculate_isotropic_distribution(u);

        // define further variables
        bool chol_flag = false;
        static thread_local VectorType g_k(dimRange, size_t(0));
        static thread_local VectorType beta_in(dimRange, size_t(0));
        static thread_local VectorType beta_out(dimRange, size_t(0));
        static thread_local MatrixType T_k(dimRange, dimRange, size_t(0));

        const auto r_max = r_sequence_.back();
        for (const auto& r : r_sequence_) {
          beta_in.deep_copy((beta_cache_ && !boundary) ? *beta_cache_ : ((beta_cache_boundary_ && boundary)
                                                                             ? *beta_cache_boundary_
                                                                             : VectorType(alpha_iso, true, size_t(0))));
          T_k.deep_copy((T_cache_ && !boundary) ? *T_cache_ : ((T_cache_boundary_ && boundary) ? *T_cache_boundary_
                                                                                               : T_minus_one_));
          // normalize u
          VectorType r_times_u_iso(u_iso, true, size_t(0));
          r_times_u_iso *= r;
          VectorType v(u, true, size_t(0));
          v *= 1 - r;
          v += r_times_u_iso;
          // calculate T_k u
          static thread_local VectorType v_k(dimRange, size_t(0));
          XT::LA::solve_lower_triangular(T_k, v_k, v);
          // calculate values of basis p = S_k m
          static thread_local BasisValuesMatrixType P_k(M_.size(), VectorType(dimRange, size_t(0)));
          for (size_t ii = 0; ii < M_.size(); ++ii)
            XT::LA::solve_lower_triangular(T_k, P_k[ii], M_[ii]);
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
            static thread_local StateRangeType error(0);
            static thread_local VectorType tmp_m(dimRange, size_t(0));
            static thread_local VectorType Tinv_m(dimRange, size_t(0));
            const auto& tmp_m_entries = tmp_m.entries();
            const auto& tmp_m_indices = tmp_m.indices();
            std::fill(error.begin(), error.end(), 0.);
            for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
              tmp_m.deep_copy(M_[ll]);
              XT::LA::solve_lower_triangular(T_k, Tinv_m, tmp_m);
              tmp_m *= std::exp(beta_out * Tinv_m) * quadrature_[ll].weight();
              for (size_t nn = 0; nn < tmp_m_indices.size(); ++nn)
                error[tmp_m_indices[nn]] += tmp_m_entries[nn];
            }
            for (size_t nn = 0; nn < v.indices().size(); ++nn)
              error[v.indices()[nn]] -= v.entries()[nn];
            // calculate descent direction d_k;
            static thread_local VectorType d_k = g_k;
            d_k.deep_copy(g_k);
            d_k *= -1;
            static thread_local VectorType T_k_inv_transp_d_k(dimRange, size_t(0));
            try {
              XT::LA::solve_lower_triangular_transposed(T_k, T_k_inv_transp_d_k, d_k);
            } catch (const Dune::FMatrixError& e) {
              if (r < r_max)
                break;
              else
                DUNE_THROW(Dune::FMatrixError, e.what());
            }
            if (error.two_norm() < tau_ && std::exp(5 * T_k_inv_transp_d_k.l1_norm()) < 1 + epsilon_gamma_) {
              XT::LA::solve_lower_triangular_transposed(T_k, alpha, beta_out);
              goto outside_all_loops;
            } else {
              RangeFieldType zeta_k = 1;
              beta_in.deep_copy(beta_out);
              // backtracking line search
              while (zeta_k > epsilon_ * beta_out.l2_norm() / d_k.l2_norm()) {
                RangeFieldType f(0);
                static thread_local VectorType beta_new = d_k;
                beta_new.deep_copy(d_k);
                beta_new *= zeta_k;
                beta_new += beta_out;
                for (size_t ll = 0; ll < quadrature_.size(); ++ll)
                  f += quadrature_[ll].weight() * std::exp(beta_new * P_k[ll]);
                f -= beta_new * v_k;
                if (XT::Common::FloatCmp::le(f, f_k + xi_ * zeta_k * (g_k * d_k))) {
                  beta_in.deep_copy(beta_new);
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
          if (!alpha_cache_) {
            alpha_cache_ = std::make_unique<std::pair<StateRangeType, VectorType>>(std::make_pair(u, alpha));
            beta_cache_ = std::make_unique<VectorType>(beta_out);
            T_cache_ = std::make_unique<MatrixType>(T_k);
          } else {
            alpha_cache_->first = u;
            alpha_cache_->second.deep_copy(alpha);
            beta_cache_->deep_copy(beta_out);
            T_cache_->deep_copy(T_k);
          }
        } else {
          if (!alpha_cache_boundary_) {
            alpha_cache_boundary_ = std::make_unique<std::pair<StateRangeType, VectorType>>(std::make_pair(u, alpha));
            beta_cache_boundary_ = std::make_unique<VectorType>(beta_out);
            T_cache_boundary_ = std::make_unique<MatrixType>(T_k);
          } else {
            alpha_cache_boundary_->first = u;
            alpha_cache_boundary_->second.deep_copy(alpha);
            beta_cache_boundary_->deep_copy(beta_out);
            T_cache_boundary_->deep_copy(T_k);
          }
        } // else (!boundary)
        mutex_.unlock();
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
        const auto& m = M_[ll];
        const auto& m_entries = m.entries();
        const auto& m_indices = m.indices();
        const auto& position = quadrature_[ll].position();
        const auto& weight = quadrature_[ll].weight();
        const auto factor = std::exp(alpha * M_[ll]) * weight * position[col];
        for (size_t kk = 0; kk < m_entries.size(); ++kk)
          ret[m_indices[kk]] += m_entries[kk] * factor;
      } // ll
    } // void evaluate_col(...)


    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      thread_local FieldMatrixType H;
      calculate_hessian(alpha, M_, H);
      helper<dimDomain>::partial_u(alpha, M_, H, ret, this);
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      const auto alpha = get_alpha(x_local, u, param);
      thread_local FieldMatrixType H;
      calculate_hessian(alpha, M_, H);
      helper<dimDomain>::partial_u_col(col, alpha, M_, H, ret, this);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    template <size_t domainDim = dimDomain, class anything = void>
    struct helper
    {
      static void axpy(const size_t dd, RangeType& ret, const RangeFieldType& factor, const VectorType& m)
      {
        const auto& m_entries = m.entries();
        const auto& m_indices = m.indices();
        for (size_t kk = 0; kk < m_entries.size(); ++kk)
          ret[m_indices[kk]][dd] += factor * m_entries[kk];
      }

      static void partial_u(const VectorType& alpha,
                            const BasisValuesMatrixType& M,
                            const FieldMatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        for (size_t dd = 0; dd < domainDim; ++dd) {
          entropy_flux->calculate_J(alpha, M, ret[dd], dd);
          calculate_A_Binv(ret[dd], H, dd > 0);
        }
      } // void partial_u(...)

      static void partial_u_col(const size_t col,
                                const VectorType& alpha,
                                const BasisValuesMatrixType& M,
                                const FieldMatrixType& H,
                                ColPartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(alpha, M, ret, col);
        calculate_A_Binv(ret, H);
      } // void partial_u_col(...)
    }; // class helper<...>

    template <class anything>
    struct helper<1, anything>
    {
      static void axpy(const size_t dd, RangeType& ret, const RangeFieldType& factor, const VectorType& m)
      {
        assert(dd == 0);
        const auto& m_entries = m.entries();
        const auto& m_indices = m.indices();
        for (size_t kk = 0; kk < m_entries.size(); ++kk)
          ret[m_indices[kk]] += factor * m_entries[kk];
      }

      static void partial_u(const VectorType& alpha,
                            const BasisValuesMatrixType& M,
                            FieldMatrixType& H,
                            PartialURangeType& ret,
                            const Localfunction* entropy_flux)
      {
        entropy_flux->calculate_J(alpha, M, ret, 0);
        calculate_A_Binv(ret, H);
      } // void partial_u(...)


      static void partial_u_col(const size_t col,
                                const VectorType& alpha,
                                const BasisValuesMatrixType& M,
                                FieldMatrixType& H,
                                PartialURangeType& ret,
                                const Localfunction* entropy_flux)
      {
        assert(col == 0);
        partial_u(alpha, M, H, ret, entropy_flux);
      } // void partial_u(...)
    }; // class helper<1, ...>

    // calculates A = A B^{-1}. B is assumed to be symmetric positive definite.
    static void calculate_A_Binv(FieldMatrixType& A, const FieldMatrixType& B, bool reuse_L = false)
    {
      // if B = LL^T, then we have to calculate ret = A (L^T)^{-1} L^{-1} = C L^{-1}
      static thread_local MatrixType L(dimRange, dimRange, size_t(0));
      if (!reuse_L) {
        // calculate B = LL^T
        bool chol_flag = XT::LA::cholesky_L(B, L);
        if (!chol_flag)
          DUNE_THROW(Dune::MathError, "B has to be symmetric positive definite!");
      }
      thread_local VectorType tmp_vec;
      for (size_t ii = 0; ii < dimRange; ++ii) {
        // calculate C = A (L^T)^{-1} and store in B
        XT::LA::solve_lower_triangular(L, tmp_vec, A[ii]);
        // calculate ret = C L^{-1}
        XT::LA::solve_lower_triangular_transposed(L, A[ii], tmp_vec);
      } // ii
    } // void calculate_A_Binv(...)

    void calculate_hessian(const VectorType& alpha, const BasisValuesMatrixType& M, FieldMatrixType& H) const
    {
      std::fill(H.begin(), H.end(), 0);
      static thread_local VectorType m_times_factor(dimRange, size_t(0));
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        m_times_factor.deep_copy(M[ll]);
        m_times_factor *= std::exp(alpha * M[ll]) * quadrature_[ll].weight();
        const auto& m_times_factor_entries = m_times_factor.entries();
        const auto& m_entries = M[ll].entries();
        const auto& m_indices = M[ll].indices();
        // add m m^T * factor
        for (size_t kk = 0; kk < m_entries.size(); ++kk)
          for (size_t nn = 0; nn < m_entries.size(); ++nn) {
            // matrix is symmetric, we only use lower triangular part
            if (m_indices[nn] > m_indices[kk])
              break;
            H[m_indices[kk]][m_indices[nn]] += m_entries[kk] * m_times_factor_entries[nn];
          }
      } // quadrature points for loop
    } // void calculate_hessian(...)

    // J = df/dalpha is the derivative of the flux with respect to alpha.
    // As F = (f_1, f_2, f_3) is matrix-valued
    // (div f = \sum_{i=1}^d \partial_{x_i} f_i  = \sum_{i=1}^d \partial_{x_i} < v_i m \hat{psi}(alpha) > is
    // vector-valued),
    // the derivative is the vector of matrices (df_1/dalpha, df_2/dalpha, ...)
    // this function returns the dd-th matrix df_dd/dalpha of J
    void calculate_J(const VectorType& alpha,
                     const BasisValuesMatrixType& M,
                     Dune::FieldMatrix<double, dimRange, StateType::dimRange>& J_dd,
                     const size_t dd) const
    {
      assert(dd < dimRangeCols);
      std::fill(J_dd.begin(), J_dd.end(), 0);
      static thread_local VectorType m_times_factor(dimRange, size_t(0));
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        const auto& v = quadrature_[ll].position();
        m_times_factor.deep_copy(M[ll]);
        m_times_factor *= v[dd] * std::exp(alpha * M[ll]) * quadrature_[ll].weight();
        const auto& m_times_factor_entries = m_times_factor.entries();
        const auto& m_times_factor_indices = m_times_factor.indices();
        const auto& m_entries = M[ll].entries();
        const auto& m_indices = M[ll].indices();
        // add m m^T * factor
        for (size_t kk = 0; kk < m_entries.size(); ++kk)
          for (size_t nn = 0; nn < m_times_factor_entries.size(); ++nn)
            J_dd[m_indices[kk]][m_times_factor_indices[nn]] += m_entries[kk] * m_times_factor_entries[nn];
      } // quadrature points for loop
      //       symmetric update for lower triangular part of J
      //      for (size_t rr = 0; rr < dimRange; ++rr)
      //        for (size_t cc = 0; cc < rr; ++cc)
      //          J_dd[rr][cc] = J_dd[cc][rr];
    } // void calculate_J(...)

    void change_basis(bool& chol_flag,
                      const VectorType& beta_in,
                      VectorType& v_k,
                      BasisValuesMatrixType& P_k,
                      MatrixType& T_k,
                      VectorType& g_k,
                      VectorType& beta_out) const
    {
      thread_local auto H = XT::Common::make_unique<FieldMatrixType>(0);
      thread_local MatrixType L(dimRange, dimRange, size_t(0));

      calculate_hessian(beta_in, P_k, *H);
      chol_flag = XT::LA::cholesky_L(*H, L);
      if (chol_flag == false)
        return;
      static thread_local VectorType Linv_P(dimRange, size_t(0));
      for (size_t ll = 0; ll < P_k.size(); ++ll) {
        XT::LA::solve_lower_triangular(L, Linv_P, P_k[ll]);
        P_k[ll].deep_copy(Linv_P);
      }
      T_k.rightmultiply(L);
      L.mtv(beta_in, beta_out);

      VectorType& v_k_tmp = Linv_P;
      XT::LA::solve_lower_triangular(L, v_k_tmp, v_k);
      v_k.deep_copy(v_k_tmp);
      StateRangeType g_k_tmp(0);
      for (size_t kk = 0; kk < v_k.indices().size(); ++kk)
        g_k_tmp[v_k.indices()[kk]] = -v_k.entries()[kk];
      VectorType& contribution = Linv_P;
      for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
        // assumes that the first basis function is constant
        //      g_k[0] += std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight() * P_k[0][0];
        contribution.deep_copy(P_k[ll]);
        contribution *= std::exp(beta_out * P_k[ll]) * quadrature_[ll].weight();
        g_k_tmp += contribution;
      }
      g_k = g_k_tmp;
    } // void change_basis(...)

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
    std::unique_ptr<std::pair<StateRangeType, VectorType>>& alpha_cache_;
    std::unique_ptr<VectorType>& beta_cache_;
    std::unique_ptr<MatrixType>& T_cache_;
    std::unique_ptr<std::pair<StateRangeType, VectorType>>& alpha_cache_boundary_;
    std::unique_ptr<VectorType>& beta_cache_boundary_;
    std::unique_ptr<MatrixType>& T_cache_boundary_;
    std::mutex& mutex_;
  }; // class Localfunction>

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
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
                                           T_cache_[index_set_.size(0) + index],
                                           mutexes_[index]);
  }


  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param);
    const auto alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor);
    StateRangeType ret(0);
    for (size_t ll = 0; ll < quadrature_.size(); ++ll) {
      const auto& position = quadrature_[ll].position();
      const auto& weight = quadrature_[ll].weight();
      const auto& m = M_[ll];
      const auto factor = position * n_ij > 0 ? std::exp(alpha_i * m) * weight : std::exp(alpha_j * m) * weight;
      auto contribution = m;
      contribution *= position[dd] * factor * n_ij[dd];
      ret += contribution;
    }
    return ret;
  } // StateRangeType evaluate_kinetic_flux(...)

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
  const MatrixType& T_minus_one_;
  const std::string name_;
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see
  // constructor)
  mutable std::vector<std::unique_ptr<std::pair<StateRangeType, VectorType>>> alpha_cache_;
  mutable std::vector<std::unique_ptr<VectorType>> beta_cache_;
  mutable std::vector<std::unique_ptr<MatrixType>> T_cache_;
  mutable std::vector<std::mutex> mutexes_;
};

/**
 * Specialization of EntropyBasedLocalFlux for 1D Hatfunctions
 */
template <class GridLayerImp, class U>
class EntropyBasedLocalFlux<Hyperbolic::Problems::HatFunctions<typename U::DomainFieldType,
                                                               1,
                                                               typename U::RangeFieldType,
                                                               U::dimRange,
                                                               1,
                                                               1>,
                            GridLayerImp,
                            U,
                            1>
    : public XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                             typename U::DomainFieldType,
                                                             GridLayerImp::dimension,
                                                             U,
                                                             0,
                                                             typename U::RangeFieldType,
                                                             U::dimRange,
                                                             1>
{
  typedef typename XT::Functions::LocalizableFluxFunctionInterface<typename GridLayerImp::template Codim<0>::Entity,
                                                                   typename U::DomainFieldType,
                                                                   GridLayerImp::dimension,
                                                                   U,
                                                                   0,
                                                                   typename U::RangeFieldType,
                                                                   U::dimRange,
                                                                   1>
      BaseType;
  typedef EntropyBasedLocalFlux ThisType;

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
  typedef Hyperbolic::Problems::HatFunctions<DomainFieldType, 1, RangeFieldType, dimRange, 1, 1> BasisfunctionType;
  typedef GridLayerImp GridLayerType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> QuadratureRuleType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

  explicit EntropyBasedLocalFlux(
      const BasisfunctionType& basis_functions,
      const GridLayerType& grid_layer,
      const QuadratureRuleType& /*quadrature*/,
      const RangeFieldType tau = 1e-9,
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
    : index_set_(grid_layer.indexSet())
    , v_points_(basis_functions.triangulation())
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
    , alpha_cache_(2 * index_set_.size(0))
    , mutexes_(index_set_.size(0))
  {
  }

  class Localfunction : public LocalfunctionType
  {
  public:
    using LocalfunctionType::dimDomain;
    using typename LocalfunctionType::ColRangeType;
    using typename LocalfunctionType::ColPartialURangeType;

    Localfunction(const EntityType& e,
                  const RangeType& v_points,
                  const RangeFieldType tau,
                  const RangeFieldType epsilon_gamma,
                  const RangeFieldType chi,
                  const RangeFieldType xi,
                  const std::vector<RangeFieldType>& r_sequence,
                  const size_t k_0,
                  const size_t k_max,
                  const RangeFieldType epsilon,
                  const RangeFieldType taylor_tol,
                  const size_t taylor_order,
                  std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache,
                  std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_boundary,
                  std::mutex& mutex)
      : LocalfunctionType(e)
      , v_points_(v_points)
      , tau_(tau)
      , epsilon_gamma_(epsilon_gamma)
      , chi_(chi)
      , xi_(xi)
      , r_sequence_(r_sequence)
      , k_0_(k_0)
      , k_max_(k_max)
      , epsilon_(epsilon)
      , taylor_order_(taylor_order)
      , taylor_tol_(taylor_tol)
      , alpha_cache_(alpha_cache)
      , alpha_cache_boundary_(alpha_cache_boundary)
      , mutex_(mutex)
    {
    }

    using LocalfunctionType::entity;

    RangeType get_alpha(const DomainType& /*x_local*/, const StateRangeType& u, const XT::Common::Parameter param) const
    {
      const bool boundary = bool(param.get("boundary")[0]);
      RangeType alpha(0.);

      mutex_.lock();
      if (!boundary && alpha_cache_ && XT::Common::FloatCmp::eq(alpha_cache_->first, u)) {
        alpha = alpha_cache_->second;
        mutex_.unlock();
        return alpha;
      } else if (boundary && alpha_cache_boundary_ && XT::Common::FloatCmp::eq(alpha_cache_boundary_->first, u)) {
        alpha = alpha_cache_boundary_->second;
        mutex_.unlock();
        return alpha;
      } else {
        // if value has already been calculated for this entity at this time, skip computation

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
          RangeType alpha_k = (alpha_cache_ && !boundary)
                                  ? alpha_cache_->second
                                  : ((alpha_cache_boundary_ && boundary) ? alpha_cache_boundary_->second : alpha_iso);
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
                      (v_points_[nn] - v_points_[nn - 1]) * ((std::exp(alpha_k[nn]) + std::exp(alpha_k[nn - 1]))
                                                                 / std::pow(alpha_k[nn] - alpha_k[nn - 1], 2)
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
                      (v_points_[nn + 1] - v_points_[nn]) * ((std::exp(alpha_k[nn + 1]) + std::exp(alpha_k[nn]))
                                                                 / std::pow(alpha_k[nn + 1] - alpha_k[nn], 2)
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
        if (!boundary) {
          if (!alpha_cache_) {
            alpha_cache_ = std::make_unique<std::pair<StateRangeType, StateRangeType>>(std::make_pair(u, alpha));
          } else {
            alpha_cache_->first = u;
            alpha_cache_->second = alpha;
          }
        } else {
          if (!alpha_cache_boundary_) {
            alpha_cache_boundary_ =
                std::make_unique<std::pair<StateRangeType, StateRangeType>>(std::make_pair(u, alpha));
          } else {
            alpha_cache_boundary_->first = u;
            alpha_cache_boundary_->second = alpha;
          }
        } // else (!boundary)
        mutex_.unlock();
      } // else ( value has not been calculated before )

      return alpha;
    } // ... get_alpha(...)

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

      // calculate < \mu m G_\alpha(u) >
      ret *= 0.;
      for (size_t nn = 0; nn < dimRange; ++nn) {
        if (nn > 0) {
          if (XT::Common::FloatCmp::ne(alpha[nn], alpha[nn - 1], taylor_tol_)) {
            ret[nn] += 2. * std::pow(v_points_[nn] - v_points_[nn - 1], 2) / std::pow(alpha[nn] - alpha[nn - 1], 3)
                           * (std::exp(alpha[nn]) - std::exp(alpha[nn - 1]))
                       + (v_points_[nn] - v_points_[nn - 1]) / std::pow(alpha[nn] - alpha[nn - 1], 2)
                             * (v_points_[nn - 1] * (std::exp(alpha[nn]) + std::exp(alpha[nn - 1]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       + v_points_[nn] * (v_points_[nn] - v_points_[nn - 1]) / (alpha[nn] - alpha[nn - 1])
                             * std::exp(alpha[nn]);
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
            ret[nn] += -2. * std::pow(v_points_[nn + 1] - v_points_[nn], 2) / std::pow(alpha[nn + 1] - alpha[nn], 3)
                           * (std::exp(alpha[nn + 1]) - std::exp(alpha[nn]))
                       + (v_points_[nn + 1] - v_points_[nn]) / std::pow(alpha[nn + 1] - alpha[nn], 2)
                             * (v_points_[nn + 1] * (std::exp(alpha[nn + 1]) + std::exp(alpha[nn]))
                                - 2 * v_points_[nn] * std::exp(alpha[nn]))
                       - v_points_[nn] * (v_points_[nn + 1] - v_points_[nn]) / (alpha[nn + 1] - alpha[nn])
                             * std::exp(alpha[nn]);
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
    } // void evaluate(...)

    virtual void evaluate_col(const size_t col,
                              const DomainType& x_local,
                              const StateRangeType& u,
                              ColRangeType& ret,
                              const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      evaluate(x_local, u, ret, param);
    } // void evaluate_col(...)


    virtual void partial_u(const DomainType& x_local,
                           const StateRangeType& u,
                           PartialURangeType& ret,
                           const XT::Common::Parameter& param) const override
    {
      DUNE_THROW(NotImplemented, "");
    }

    virtual void partial_u_col(const size_t col,
                               const DomainType& x_local,
                               const StateRangeType& u,
                               ColPartialURangeType& ret,
                               const XT::Common::Parameter& param) const override
    {
      assert(col == 0);
      partial_u(x_local, u, ret, param);
    }

    static std::string static_id()
    {
      return "gdt.entropybasedlocalflux";
    }

  private:
    const RangeType& v_points_;
    const RangeFieldType tau_;
    const RangeFieldType epsilon_gamma_;
    const RangeFieldType chi_;
    const RangeFieldType xi_;
    const std::vector<RangeFieldType>& r_sequence_;
    const size_t k_0_;
    const size_t k_max_;
    const RangeFieldType epsilon_;
    const RangeFieldType taylor_tol_;
    const size_t taylor_order_;
    const std::string name_;
    std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_;
    std::unique_ptr<std::pair<StateRangeType, StateRangeType>>& alpha_cache_boundary_;
    std::mutex& mutex_;
  }; // class Localfunction>

  static std::string static_id()
  {
    return "gdt.entropybasedflux";
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return derived_local_function(entity);
  }

  std::unique_ptr<Localfunction> derived_local_function(const EntityType& entity) const
  {
    const auto& index = index_set_.index(entity);
    return std::make_unique<Localfunction>(entity,
                                           v_points_,
                                           tau_,
                                           epsilon_gamma_,
                                           chi_,
                                           xi_,
                                           r_sequence_,
                                           k_0_,
                                           k_max_,
                                           epsilon_,
                                           taylor_tol_,
                                           taylor_order_,
                                           alpha_cache_[index],
                                           alpha_cache_[index_set_.size(0) + index],
                                           mutexes_[index]);
  }

  // calculate \sum_{i=1}^d < v_i m \psi > n_i, where n is the unit outer normal,
  // m is the basis function vector, phi_u is the ansatz corresponding to u
  // and x, v, t are the space, velocity and time variable, respectively
  // As we are using cartesian grids, n_i == 0 in all but one dimension, so only evaluate for i == dd
  StateRangeType evaluate_kinetic_flux(const EntityType& entity,
                                       const DomainType& x_local_entity,
                                       const StateRangeType& u_i,
                                       const EntityType& neighbor,
                                       const DomainType& x_local_neighbor,
                                       const StateRangeType u_j,
                                       const DomainType& n_ij,
                                       const size_t dd,
                                       const XT::Common::Parameter& param,
                                       const XT::Common::Parameter& param_neighbor) const
  {
    assert(v_points_.size() % 2 && "Not implemented for odd number of points!");
    assert(dd == 0);
    // calculate < \mu m G_\alpha(u) > * n_ij
    const auto local_function_entity = derived_local_function(entity);
    const auto local_function_neighbor = derived_local_function(neighbor);
    const auto alpha_i = local_function_entity->get_alpha(x_local_entity, u_i, param);
    const auto alpha_j = local_function_neighbor->get_alpha(x_local_neighbor, u_j, param_neighbor);
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
  } // StateRangeType evaluate_kinetic_flux(...)

  const typename GridLayerType::IndexSet& index_set_;
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
  // Use unique_ptr in the vectors to avoid the memory cost for storing twice as many matrices or vectors as needed
  // (see constructor)
  mutable std::vector<std::unique_ptr<std::pair<StateRangeType, StateRangeType>>> alpha_cache_;
  mutable std::vector<std::mutex> mutexes_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
