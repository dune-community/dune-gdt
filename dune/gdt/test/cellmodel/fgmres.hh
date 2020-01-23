// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GDT_TEST_CELLMODEL_FGMRES
#define DUNE_GDT_TEST_CELLMODEL_FGMRES

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/conditional.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/timer.hh>

#include <dune/istl/allocator.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/eigenvalue/arpackpp.hh>
#include <dune/istl/istlexception.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solver.hh>

namespace Dune {

/**
   \brief implements the Generalized Minimal Residual (GMRes) method

   GMRes solves the unsymmetric linear system Ax = b using the
   Generalized Minimal Residual method as described the SIAM Templates
   book (http://www.netlib.org/templates/templates.pdf).

   \tparam X trial vector, vector type of the solution
   \tparam Y test vector, vector type of the RHS
   \tparam F vector type for orthonormal basis of Krylov space
 */

template <class X, class Y = X, class F = Y>
class FGMResSolver : public IterativeSolver<X, Y>
{
public:
  using typename IterativeSolver<X, Y>::domain_type;
  using typename IterativeSolver<X, Y>::range_type;
  using typename IterativeSolver<X, Y>::field_type;
  using typename IterativeSolver<X, Y>::real_type;

private:
  using typename IterativeSolver<X, X>::scalar_real_type;

  //! \bief field_type Allocator retrieved from domain type
  using fAlloc = ReboundAllocatorType<X, field_type>;
  //! \bief real_type Allocator retrieved from domain type
  using rAlloc = ReboundAllocatorType<X, real_type>;

public:
  /*!
     \brief Set up FGMResSolver solver.

     \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
     \param restart number of GMRes cycles before restart
   */
  FGMResSolver(LinearOperator<X, Y>& op,
               Preconditioner<X, Y>& prec,
               scalar_real_type reduction,
               int restart,
               int maxit,
               int verbose)
    : IterativeSolver<X, Y>::IterativeSolver(op, prec, reduction, maxit, verbose)
    , _restart(restart)
  {}

  /*!
     \brief Set up FGMResSolver solver.

     \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
     \param restart number of GMRes cycles before restart
   */
  FGMResSolver(LinearOperator<X, Y>& op,
               ScalarProduct<X>& sp,
               Preconditioner<X, Y>& prec,
               scalar_real_type reduction,
               int restart,
               int maxit,
               int verbose)
    : IterativeSolver<X, Y>::IterativeSolver(op, sp, prec, reduction, maxit, verbose)
    , _restart(restart)
  {}

  /*!
     \brief Apply inverse operator.

     \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)

     \note Currently, the FGMResSolver aborts when it detects a
           breakdown.
   */
  virtual void apply(X& x, Y& b, InverseOperatorResult& res)
  {
    apply(x, b, max_value(_reduction), res);
  }

  /*!
     \brief Apply inverse operator.

     \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)

     \note Currently, the FGMResSolver aborts when it detects a
           breakdown.
   */
  virtual void apply(X& x, Y& b, double reduction, InverseOperatorResult& res)
  {
    using std::abs;
    const real_type EPSILON = 1e-80;
    const int m = _restart;
    real_type norm, norm_old = 0.0, norm_0;
    int j = 1;
    std::vector<field_type, fAlloc> s(m + 1), sn(m);
    std::vector<real_type, rAlloc> cs(m);
    // need copy of rhs if GMRes has to be restarted
    Y b2(b);
    // helper vector
    Y w(b);
    std::vector<std::vector<field_type, fAlloc>> H(m + 1, s);
    std::vector<F> v(m + 1, b);
    std::vector<F> z(m + 1, b);

    // start timer
    Dune::Timer watch;
    watch.reset();

    // clear solver statistics and set res.converged to false
    res.clear();
    _prec->pre(x, b);

    // calculate defect and overwrite rhs with it
    _op->applyscaleadd(-1.0, x, b); // b -= Ax
    v[0] = b;
    norm_0 = _sp->norm(v[0]); // beta
    norm = norm_0;
    norm_old = norm;

    // print header
    if (_verbose > 0) {
      std::cout << "=== FGMResSolver" << std::endl;
      if (_verbose > 1) {
        this->printHeader(std::cout);
        this->printOutput(std::cout, 0, norm_0);
      }
    }

    if (all_true(norm_0 < EPSILON)) {
      _prec->post(x);
      res.converged = true;
      if (_verbose > 0) // final print
        print_result(res);
    }

    while (j <= _maxit && res.converged != true) {

      int i = 0;
      v[0] *= real_type(1.0) / norm;
      s[0] = norm;
      for (i = 1; i < m + 1; i++)
        s[i] = 0.0;

      for (i = 0; i < m && j <= _maxit && res.converged != true; i++, j++) {
        w = 0.0;
        // do Arnoldi algorithm
        _prec->apply(z[i], v[i]);
        _op->apply(z[i], w);
        z[i + 1] = z[i]; // initial guess for preconditioner in next step
        for (int k = 0; k < i + 1; k++) {
          // notice that _sp->dot(v[k],w) = v[k]\adjoint w
          // so one has to pay attention to the order
          // in the scalar product for the complex case
          // doing the modified Gram-Schmidt algorithm
          H[k][i] = _sp->dot(v[k], w);
          // w -= H[k][i] * v[k]
          w.axpy(-H[k][i], v[k]);
        }
        H[i + 1][i] = _sp->norm(w);
        if (all_true(abs(H[i + 1][i]) < EPSILON))
          DUNE_THROW(SolverAbort, "breakdown in GMRes - |w| == 0.0 after " << j << " iterations");

        // normalize new vector
        v[i + 1] = w;
        v[i + 1] *= real_type(1.0) / H[i + 1][i];

        // update QR factorization
        for (int k = 0; k < i; k++)
          applyPlaneRotation(H[k][i], H[k + 1][i], cs[k], sn[k]);

        // compute new givens rotation
        generatePlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
        // finish updating QR factorization
        applyPlaneRotation(H[i][i], H[i + 1][i], cs[i], sn[i]);
        applyPlaneRotation(s[i], s[i + 1], cs[i], sn[i]);

        // norm of the defect is the last component the vector s
        norm = abs(s[i + 1]);

        // print current iteration statistics
        if (_verbose > 1) {
          this->printOutput(std::cout, j, norm, norm_old);
        }

        norm_old = norm;

        // check convergence
        if (all_true(norm < real_type(reduction) * norm_0))
          res.converged = true;

      } // end for

      // calculate update vector
      w = 0.0;
      update(w, i, H, s, z);
      // and current iterate
      x += w;

      // restart GMRes if convergence was not achieved,
      // i.e. linear defect has not reached desired reduction
      // and if j < _maxit
      if (res.converged != true && j <= _maxit) {

        if (_verbose > 0)
          std::cout << "=== GMRes::restart" << std::endl;
        // get saved rhs
        b = b2;
        // calculate new defect
        _op->applyscaleadd(-1.0, x, b); // b -= Ax;
        // calculate preconditioned defect
        v[0] = b;
        norm = _sp->norm(v[0]);
        norm_old = norm;
      }

    } // end while

    // postprocess preconditioner
    _prec->post(x);

    // save solver statistics
    res.iterations = j - 1; // it has to be j-1!!!
    res.reduction = static_cast<double>(max_value(norm / norm_0));
    res.conv_rate = pow(res.reduction, 1.0 / (j - 1));
    res.elapsed = watch.elapsed();

    if (_verbose > 0)
      print_result(res);
  }

private:
  void print_result(const InverseOperatorResult& res) const
  {
    int k = res.iterations > 0 ? res.iterations : 1;
    std::cout << "=== rate=" << res.conv_rate << ", T=" << res.elapsed << ", TIT=" << res.elapsed / k
              << ", IT=" << res.iterations << std::endl;
  }

  void update(X& w,
              int i,
              const std::vector<std::vector<field_type, fAlloc>>& H,
              const std::vector<field_type, fAlloc>& s,
              const std::vector<X>& z)
  {
    // solution vector of the upper triangular system
    std::vector<field_type, fAlloc> y(s);

    // backsolve
    for (int a = i - 1; a >= 0; a--) {
      field_type rhs(s[a]);
      for (int b = a + 1; b < i; b++)
        rhs -= H[a][b] * y[b];
      y[a] = rhs / H[a][a];

      // compute update on the fly
      // w += y[a]*z[a]
      w.axpy(y[a], z[a]);
    }
  }

  template <typename T>
  typename std::enable_if<std::is_same<field_type, real_type>::value, T>::type conjugate(const T& t)
  {
    return t;
  }

  template <typename T>
  typename std::enable_if<!std::is_same<field_type, real_type>::value, T>::type conjugate(const T& t)
  {
    using std::conj;
    return conj(t);
  }

  // helper function to extract the real value of a real or complex number
  inline real_type to_real(const real_type& v)
  {
    return v;
  }

  inline real_type to_real(const std::complex<real_type>& v)
  {
    return v.real();
  }

  void generatePlaneRotation(field_type& dx, field_type& dy, real_type& cs, field_type& sn)
  {
    using std::abs;
    using std::max;
    using std::min;
    using std::sqrt;
    const real_type eps = 1e-15;
    real_type norm_dx = abs(dx);
    real_type norm_dy = abs(dy);
    real_type norm_max = max(norm_dx, norm_dy);
    real_type norm_min = min(norm_dx, norm_dy);
    real_type temp = norm_min / norm_max;
    // we rewrite the code in a vectorizable fashion
    cs = cond(norm_dy < eps,
              real_type(1.0),
              cond(norm_dx < eps,
                   real_type(0.0),
                   cond(norm_dy > norm_dx,
                        real_type(1.0) / sqrt(real_type(1.0) + temp * temp) * temp,
                        real_type(1.0) / sqrt(real_type(1.0) + temp * temp))));
    sn = cond(norm_dy < eps,
              field_type(0.0),
              cond(norm_dx < eps,
                   field_type(1.0),
                   cond(norm_dy > norm_dx,
                        field_type(1.0) / sqrt(real_type(1.0) + temp * temp) * dx * conjugate(dy) / norm_dx / norm_dy,
                        field_type(1.0) / sqrt(real_type(1.0) + temp * temp) * conjugate(dy / dx))));
  }


  void applyPlaneRotation(field_type& dx, field_type& dy, real_type& cs, field_type& sn)
  {
    field_type temp = cs * dx + sn * dy;
    dy = -conjugate(sn) * dx + cs * dy;
    dx = temp;
  }

  using IterativeSolver<X, Y>::_op;
  using IterativeSolver<X, Y>::_prec;
  using IterativeSolver<X, Y>::_sp;
  using IterativeSolver<X, Y>::_reduction;
  using IterativeSolver<X, Y>::_maxit;
  using IterativeSolver<X, Y>::_verbose;
  int _restart;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_FGMRES
