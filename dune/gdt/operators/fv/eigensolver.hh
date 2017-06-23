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

#ifndef DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
#define DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <dune/xt/common/string.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/common.hh>

namespace Dune {
namespace GDT {


template <class LocalFluxFunctionImp,
          class VectorImp = XT::LA::CommonDenseVector<typename LocalFluxFunctionImp::RangeFieldType>,
          class MatrixImp = XT::LA::CommonDenseMatrix<typename LocalFluxFunctionImp::RangeFieldType>>
class QrHouseholderEigenSolver
{
  static const size_t stateDimRange = LocalFluxFunctionImp::StateType::dimRange;
  static const size_t dimRange = LocalFluxFunctionImp::dimRange;
  static const size_t dimRangeCols = LocalFluxFunctionImp::dimRangeCols;
  static_assert(stateDimRange == dimRange,
                "Jacobian has to be a square matrix, i.e. stateDimRange and dimRange have to be equal!");
  static const size_t rows = dimRange;
  static const size_t cols = stateDimRange;
  typedef typename LocalFluxFunctionImp::RangeFieldType RangeFieldType;
  typedef typename LocalFluxFunctionImp::DomainType DomainType;
  typedef typename LocalFluxFunctionImp::StateRangeType StateRangeType;

public:
  typedef VectorImp VectorType;
  typedef MatrixImp MatrixType;
  typedef FieldVector<VectorType, dimRangeCols> EigenValuesType;
  typedef FieldVector<MatrixType, dimRangeCols> EigenVectorsType;

public:
  QrHouseholderEigenSolver(const LocalFluxFunctionImp& local_flux_function,
                           const DomainType& x_local,
                           const StateRangeType& u,
                           const XT::Common::Parameter& param)
    : eigenvalues_(VectorType(rows))
    , eigenvectors_(MatrixType(rows, cols))
  {
    const auto partial_u = XT::Functions::JacobianRangeTypeConverter<stateDimRange, dimRange, dimRangeCols>::convert(
        local_flux_function.partial_u(x_local, u, param));


    //    MatrixType A_k, A_k_old, R_k, Q_k(rows, rows), Q(rows, rows);
    //    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
    //      RangeFieldType residual = 10.;
    //      A_k = XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange,
    //      stateDimRange>::convert(partial_u[ii]);
    //      for (size_t rr = 0; rr < rows; ++rr)
    //        Q.unit_row(rr);
    //      size_t count = 0;
    //      while (XT::Common::FloatCmp::ne(residual, 0.) && count < 100000) {
    //        A_k_old = A_k;
    //        QR_decomp(A_k, Q_k, R_k);
    //        // calculate A_{k+1} = R_k Q_k
    //        A_k = R_k * Q_k;
    //        Q = Q * Q_k;
    //        residual = (A_k_old - A_k).sup_norm();
    //        ++count;
    //      } // while(residual != 0)

    //      if (count > 90000) {
    //        std::cout << XT::Common::to_string(A_k_old) << std::endl;
    //        std::cout << XT::Common::to_string(R_k * Q_k) << std::endl;
    //        std::cout << XT::Common::to_string(Q_k * R_k) << std::endl;
    //        std::cout << XT::Common::to_string(Q_k) << std::endl;
    //        std::cout << "Iteration: " << count << ", residual: " << XT::Common::to_string(residual, 15) << std::endl;
    //        std::cout << "Partial u: " << XT::Common::to_string(partial_u) << std::endl;
    //      }
    //      // eigenvalues are the diagonal elements of A_{final}
    //      for (size_t rr = 0; rr < rows; ++rr)
    //        eigenvalues_[ii].set_entry(rr, A_k.get_entry(rr, rr));
    //      // eigenvectors are the columns of Q
    //      eigenvectors_[ii] = Q;
    //    }

    static constexpr size_t max_counts = 10000;
    MatrixType A, R_k, Q_k(rows, rows), Q(rows, rows), Q_copy(rows, rows);
    const RangeFieldType tol = 1e-15;
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      for (size_t rr = 0; rr < rows; ++rr)
        Q.unit_row(rr);
      A = XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, stateDimRange>::convert(partial_u[ii]);
      hessenberg_transformation(A);
      for (size_t jj = rows - 1; jj > 0; --jj) {
        size_t num_rows = jj + 1;
        size_t num_cols = num_rows;
        RangeFieldType residual = std::abs(A.get_entry(jj, jj - 1));
        size_t kk = 0;
        while (XT::Common::FloatCmp::gt(residual,
                                        tol * (std::abs(A.get_entry(jj, jj)) + std::abs(A.get_entry(jj - 1, jj - 1))))
               && kk < max_counts) {
          auto sigma = A.get_entry(jj, jj);
          for (size_t rr = 0; rr < num_rows; ++rr)
            A.add_to_entry(rr, rr, -sigma);

          // calculate QR decomp. If jj < rows-1, A has the form [A1 A2; 0 A3],
          // so we are only calculating the QR decomposition Q1*R1 of A1.
          // Then Q = [Q1 0; 0 I], R = [R1 Q1^T*A2; 0 A3] is a QR decomp. of A.
          QR_decomp(A, Q_k, R_k, num_rows, num_cols);
          std::cout << XT::Common::to_string(A - Q_k * R_k) << std::endl;
          std::cout << XT::Common::to_string(A) << std::endl;

          // calculate A_{k+1} = R_k Q_k + sigma I. We are only interested in the diagonal
          // elements of A, and we do not reuse the other parts of A, so we are only
          // updating the upper left part.
          for (size_t rr = 0; rr < num_rows; ++rr) {
            for (size_t cc = 0; cc < num_cols; ++cc) {
              A.set_entry(rr, cc, 0.);
              for (size_t ll = 0; ll < num_rows; ++ll)
                A.add_to_entry(rr, cc, R_k.get_entry(rr, ll) * Q_k.get_entry(ll, cc));
            } // cc
          } // rr
          for (size_t rr = 0; rr < num_rows; ++rr)
            A.add_to_entry(rr, rr, sigma);
          // calculate Q = Q * Q_k. As Q_k = (Q_k' 0; 0 I) (see above), if Q = (Q1 Q2; Q3 Q4)
          // we need to calculate Q = (Q1*Q_k Q2; Q3*Q_k Q_4)
          Q_copy = Q;
          for (size_t rr = 0; rr < rows; ++rr) {
            for (size_t cc = 0; cc < num_cols; ++cc) {
              Q.set_entry(rr, cc, 0.);
              for (size_t ll = 0; ll < num_rows; ++ll)
                Q.add_to_entry(rr, cc, Q_copy.get_entry(rr, ll) * Q_k.get_entry(ll, cc));
            } // cc
          } // rr
          ++kk;
          residual = std::abs(A.get_entry(jj, jj - 1));
        } // while(residual != 0)
        if (kk >= max_counts)
          std::cerr << "Warning: Eigen solver did not converge (stopped after 10000 iterations), result may be wrong!"
                    << std::endl;
      } // jj

      // eigenvalues are the diagonal elements of A_{final}
      for (size_t rr = 0; rr < rows; ++rr)
        eigenvalues_[ii].set_entry(rr, A.get_entry(rr, rr));
      // eigenvectors are the columns of Q
      eigenvectors_[ii] = Q;
      std::cout << XT::Common::to_string(partial_u[ii]) << std::endl;
      std::cout << XT::Common::to_string(eigenvectors_[ii]) << std::endl;
      std::cout << XT::Common::to_string(eigenvalues_[ii]) << std::endl;
    } // ii
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType eigenvectors_inverse() const
  {
    // eigenvectors[ii] is an orthogonal matrix, so the inverse is just the the transposed
    EigenVectorsType ret;
    for (size_t ii = 0; ii < dimRangeCols; ++ii)
      ret[ii] = eigenvectors_[ii].transposed();
    return ret;
  }

private:
  //! \brief modified sign function returning 1 instead of 0 if the value is 0
  RangeFieldType xi(RangeFieldType val) const
  {
    return XT::Common::FloatCmp::eq(val, 0.) ? 1. : val / std::abs(val);
  }

  RangeFieldType get_norm_x(const MatrixType& A, const size_t col_index, size_t num_rows = rows)
  {
    RangeFieldType norm(0);
    for (size_t rr = col_index; rr < num_rows; ++rr)
      norm += std::pow(A.get_entry(rr, col_index), 2);
    return std::sqrt(norm);
  }

  /** \brief This is a simple QR scheme using Householder reflections.
  * The householder matrix is written as H = I - 2 v v^T, where v = u/||u|| and u = x - s ||x|| e_1, s = +-1 has the
  * opposite sign of u_1 and x is the current column of A. The matrix H is rewritten as
  * H = I - tau w w^T, where w=u/u_1 and tau = -s u_1/||x||.
  * The num_rows and num_cols parameter is used to restrict the rows and columns, the Q and R will only contain
  * the QR decomposition of A[0:num_rows, 0:num_cols].
  * \see https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections.
  * \see http://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
  * \todo The scheme does not take into account that A is in Hessenberg form, so there could be a significant
  * speedup if an according QR scheme is used (e.g. givens rotations, see https://lp.uni-goettingen.de/get/text/2138).
  */
  void QR_decomp(const MatrixType& A, MatrixType& Q, MatrixType& R, size_t num_rows = rows, size_t num_cols = cols)
  {
    assert(num_rows >= num_cols && "Not implemented yet for rows < cols!");
    assert(A.rows() >= num_rows && A.cols() >= num_cols && "A is too small!");
    assert(Q.rows() >= num_rows && Q.cols() >= num_cols && "Q is too small!");
    R = A;
    for (size_t rr = 0; rr < num_rows; ++rr)
      Q.unit_row(rr);

    VectorType w(num_rows), tau(num_rows);
    MatrixType R_copy, Q_copy;
    for (size_t jj = 0; jj < std::min(num_rows - 1, num_cols); ++jj) {
      const auto norm_x = get_norm_x(R, jj, num_rows);
      if (XT::Common::FloatCmp::gt(norm_x, 0.)) {
        const auto s = -sign(R.get_entry(jj, jj));
        const RangeFieldType u1 = R.get_entry(jj, jj) - s * norm_x;
        w.set_entry(jj, 1.);
        for (size_t rr = jj + 1; rr < num_rows; ++rr)
          w.set_entry(rr, R.get_entry(rr, jj) / u1);
        tau.set_entry(jj, -s * u1 / norm_x);

        // TODO: improve calculation, do not use O(N^3) multiplication but O(N^2) using structure of Q
        // calculate R_new = Q_k R_old and Q_new = Q_old Q_k
        R_copy = R;
        for (size_t rr = jj; rr < num_rows; ++rr)
          for (size_t cc = 0; cc < num_cols; ++cc) {
            for (size_t kk = jj; kk < num_rows; ++kk)
              R.add_to_entry(
                  rr, cc, -1. * R_copy.get_entry(kk, cc) * tau.get_entry(jj) * w.get_entry(rr) * w.get_entry(kk));
          }

        Q_copy = Q;
        for (size_t rr = 0; rr < num_rows; ++rr)
          for (size_t cc = jj; cc < num_cols; ++cc)
            for (size_t kk = jj; kk < num_cols; ++kk)
              Q.add_to_entry(
                  rr, cc, -1. * Q_copy.get_entry(rr, kk) * tau.get_entry(jj) * w.get_entry(cc) * w.get_entry(kk));
      } // if (norm_x != 0)
    } // jj

    // choose Q such that largest entry of each colum is positive
    size_t row_index = 0;
    for (size_t cc = 0; cc < num_cols; ++cc) {
      RangeFieldType max = std::numeric_limits<RangeFieldType>::lowest();
      for (size_t rr = 0; rr < num_rows; ++rr) {
        if (std::abs(Q.get_entry(rr, cc)) > max) {
          max = std::abs(Q.get_entry(rr, cc));
          row_index = rr;
        }
      } // rr
      if (XT::Common::FloatCmp::lt(Q.get_entry(row_index, cc), 0.)) {
        // scal column of Q if largest entry is negative
        // scal row of R to ensure that still A = QR
        for (size_t rr = 0; rr < num_rows; ++rr) {
          Q.set_entry(rr, cc, -Q.get_entry(rr, cc));
          R.set_entry(cc, rr, -R.get_entry(cc, rr));
        } // rr
      } // if (largest entry negative)
    } // cc
  } // void QR_decomp(...)

  //! \brief Transform A to Hessenberg form
  //! \see https://lp.uni-goettingen.de/get/text/2137
  void hessenberg_transformation(MatrixType& A) const
  {
    static_assert(rows == cols, "Hessenberg transformation needs a square matrix!");
    assert(A.rows() == rows && A.cols() == cols && "A has wrong dimensions!");
    VectorType u(rows);
    VectorType tmp_vec(rows);
    for (size_t jj = 0; jj < rows - 2; ++jj) {
      RangeFieldType gamma = 0;
      for (size_t rr = jj + 1; rr < rows; ++rr)
        gamma += std::pow(std::abs(A.get_entry(rr, jj)), 2);
      gamma = std::sqrt(gamma);
      if (XT::Common::FloatCmp::ne(gamma, 0.)) {
        for (size_t rr = jj + 1; rr < rows; ++rr)
          u.set_entry(rr, A.get_entry(rr, jj));
        u.add_to_entry(jj + 1, xi(A.get_entry(jj + 1, jj)) * gamma);
        RangeFieldType beta = 1. / (gamma * (gamma + std::abs(A.get_entry(jj + 1, jj))));
        // calculate P A P with P = diag(I_j, (I_{n-j} - beta u u*))
        // calculate P A = A - (beta u) (u^T A) first
        for (size_t cc = 0; cc < cols; ++cc) {
          tmp_vec.set_entry(cc, 0.);
          for (size_t rr = jj + 1; rr < rows; ++rr)
            tmp_vec.add_to_entry(cc, u.get_entry(rr) * A.get_entry(rr, cc));
        } // tmp_vec now contains u^T A[j+1:rows,:]
        for (size_t rr = jj + 1; rr < rows; ++rr)
          for (size_t cc = 0; cc < cols; ++cc)
            A.add_to_entry(rr, cc, -beta * u.get_entry(rr) * tmp_vec.get_entry(cc));
        // now calculate (PA) P  = PA - (PA u) (beta u^T)
        for (size_t rr = 0; rr < rows; ++rr) {
          tmp_vec.set_entry(rr, 0.);
          for (size_t cc = jj + 1; cc < cols; ++cc)
            tmp_vec.add_to_entry(rr, A.get_entry(rr, cc) * u.get_entry(cc));
        } // tmp_vec now contains PA[:,j+1:cols] u
        for (size_t rr = 0; rr < rows; ++rr)
          for (size_t cc = jj + 1; cc < cols; ++cc)
            A.add_to_entry(rr, cc, -beta * tmp_vec.get_entry(rr) * u.get_entry(cc));
      } // if (gamma != 0)
    } // jj
  } // void hessenberg_transformation(...)

  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
}; // class QrHouseholderEigenSolver<...>


#if HAVE_EIGEN

template <class LocalFluxFunctionImp, bool selfadjoint = false>
class EigenEigenSolver
{
  static const size_t stateDimRange = LocalFluxFunctionImp::StateType::dimRange;
  static const size_t dimRange = LocalFluxFunctionImp::dimRange;
  static const size_t dimRangeCols = LocalFluxFunctionImp::dimRangeCols;
  static_assert(stateDimRange == dimRange,
                "Jacobian has to be a square matrix, i.e. stateDimRange and dimRange have to be equal!");
  typedef typename LocalFluxFunctionImp::RangeFieldType RangeFieldType;
  typedef typename LocalFluxFunctionImp::DomainType DomainType;
  typedef typename LocalFluxFunctionImp::StateRangeType StateRangeType;

public:
  typedef typename XT::LA::EigenDenseVector<RangeFieldType> VectorType;
  typedef typename XT::LA::EigenDenseMatrix<RangeFieldType> MatrixType;
  typedef FieldVector<VectorType, dimRangeCols> EigenValuesType;
  typedef FieldVector<MatrixType, dimRangeCols> EigenVectorsType;

private:
  typedef typename MatrixType::BackendType EigenMatrixBackendType;
  typedef typename std::conditional<selfadjoint,
                                    ::Eigen::SelfAdjointEigenSolver<EigenMatrixBackendType>,
                                    ::Eigen::EigenSolver<EigenMatrixBackendType>>::type EigenSolverType;

public:
  EigenEigenSolver(const LocalFluxFunctionImp& local_flux_function,
                   const DomainType& x_local,
                   const StateRangeType& u,
                   const XT::Common::Parameter& param)
  {
    const auto partial_u = XT::Functions::JacobianRangeTypeConverter<stateDimRange, dimRange, dimRangeCols>::convert(
        local_flux_function.partial_u(x_local, u, param));
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      const auto partial_u_eigen =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, stateDimRange>::convert(partial_u[ii]);
      EigenSolverType eigen_solver(partial_u_eigen.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto& eigenvalues_eigen = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto& eigenvectors_eigen = eigen_solver.eigenvectors(); // <- this should be an Eigen vector of std::complex
      assert(size_t(eigenvalues_eigen.size()) == dimRange);
      assert(XT::Common::FloatCmp::eq(VectorType(eigenvalues_eigen.imag()), VectorType(dimRange, 0.)));
      assert(XT::Common::FloatCmp::eq(MatrixType(eigenvectors_eigen.imag()), MatrixType(dimRange, dimRange, 0.)));
      eigenvalues_[ii] = VectorType(eigenvalues_eigen.real());
      eigenvectors_[ii] = MatrixType(eigenvectors_eigen.real());
    }
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType eigenvectors_inverse() const
  {
    EigenVectorsType ret;
    for (size_t ii = 0; ii < dimRangeCols; ++ii)
      ret[ii] = MatrixType(eigenvectors_[ii].backend().inverse());
    return ret;
  }

private:
  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
}; // class EigenEigenSolver<...>

#else // HAVE_EIGEN

template <class LocalFluxFunctionImp, bool selfadjoint = false>
class EigenEigenSolver
{
  static_assert(AlwaysFalse<LocalFluxFunctionImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN

template <class LocalFluxFunctionImp>
//#if HAVE_EIGEN
// using DefaultEigenSolver = EigenEigenSolver<LocalFluxFunctionImp>;
//#else
using DefaultEigenSolver = QrHouseholderEigenSolver<LocalFluxFunctionImp>;
//#endif

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
