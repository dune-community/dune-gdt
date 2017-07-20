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

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/common.hh>

namespace Dune {
namespace GDT {


template <class FieldType, size_t dimRange, size_t dimRangeCols>
class QrHouseholderEigenSolver
{
  static const size_t rows = dimRange;
  static const size_t cols = dimRange;

public:
  typedef FieldMatrix<FieldType, dimRange, dimRange> MatrixType;
  typedef FieldVector<FieldType, dimRange> VectorType;
  typedef FieldVector<VectorType, dimRangeCols> EigenValuesType;
  typedef FieldVector<std::shared_ptr<MatrixType>, dimRangeCols> EigenVectorsType;
  typedef FieldVector<MatrixType, dimRangeCols> InputMatricesType;

public:
  QrHouseholderEigenSolver(InputMatricesType& matrices_in, bool calculate_eigenvectors = false)
    : calculate_eigenvectors_(calculate_eigenvectors)
  {
    initialize(matrices_in, calculate_eigenvectors);
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
  void initialize(InputMatricesType& matrices_in, const bool calculate_eigenvectors)
  {
    static constexpr size_t max_counts = 10000;
    const FieldType tol = 1e-15;
    auto R_k = XT::Common::make_unique<MatrixType>(0.);
    auto Q_k = XT::Common::make_unique<MatrixType>(0.);
    auto Q = XT::Common::make_unique<MatrixType>(0.);
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      eigenvectors_[ii] = std::make_shared<MatrixType>();
      eigenvectors_inverse_[ii] = std::make_shared<MatrixType>();
      auto& A = matrices_in[ii];
      hessenberg_transformation(A, *Q);
      for (size_t jj = rows - 1; jj > 0; --jj) {
        size_t num_rows = jj + 1;
        size_t num_cols = num_rows;
        FieldType residual = std::abs(A[jj][jj - 1]);
        size_t kk = 0;
        while (XT::Common::FloatCmp::gt(residual, tol * (std::abs(A[jj][jj]) + std::abs(A[jj - 1][jj - 1])))
               && kk < max_counts) {
          auto sigma = A[jj][jj];
          for (size_t rr = 0; rr < num_rows; ++rr)
            A[rr][rr] -= sigma;

          // calculate QR decomp. If jj < rows-1, A has the form [A1 A2; 0 A3],
          // so we are only calculating the QR decomposition Q1*R1 of A1.
          // Then Q = [Q1 0; 0 I], R = [R1 Q1^T*A2; 0 A3] is a QR decomp. of A.
          QR_decomp(A, *Q_k, *R_k, num_rows, num_cols);

          // calculate A_{k+1} = R_k Q_k + sigma I. We are only interested in the diagonal
          // elements of A, and we do not reuse the other parts of A, so we are only
          // updating the upper left part.
          for (size_t rr = 0; rr < num_rows; ++rr) {
            for (size_t cc = 0; cc < num_cols; ++cc) {
              A[rr][cc] = 0.;
              for (size_t ll = 0; ll < num_rows; ++ll)
                A[rr][cc] += R_k[rr][ll] * Q_k[ll][cc];
            } // cc
          } // rr

          // we do not need R_k anymore in this step, so use it as temporary storage
          auto& A_copy = *R_k;
          A_copy = A;
          // update upper right part
          for (size_t rr = 0; rr < num_rows; ++rr) {
            for (size_t cc = num_cols; cc < cols; ++cc) {
              A[rr][cc] = 0.;
              for (size_t ll = 0; ll < num_rows; ++ll)
                A[rr][cc] += (*Q_k)[ll][rr] * A_copy[ll][cc];
            } // cc
          } // rr

          for (size_t rr = 0; rr < num_rows; ++rr)
            A[rr][rr] += sigma;
          // calculate Q = Q * Q_k. As Q_k = (Q_k' 0; 0 I) (see above), if Q = (Q1 Q2; Q3 Q4)
          // we need to calculate Q = (Q1*Q_k Q2; Q3*Q_k Q_4)
          auto& Q_copy = *R_k;
          Q_copy = *Q;
          for (size_t rr = 0; rr < rows; ++rr) {
            for (size_t cc = 0; cc < num_cols; ++cc) {
              (*Q)[rr][cc] = 0.;
              for (size_t ll = 0; ll < num_rows; ++ll)
                (*Q)[rr][cc] += Q_copy[rr][ll] * (*Q_k)[ll][cc];
            } // cc
          } // rr
          ++kk;
          residual = std::abs(A[jj][jj - 1]);
        } // while(residual != 0)
        if (kk >= max_counts)
          std::cerr << "Warning: Eigen solver did not converge (stopped after 10000 iterations), result may be wrong!"
                    << std::endl;
      } // jj

      // eigenvalues are the diagonal elements of A_{final}
      for (size_t rr = 0; rr < rows; ++rr)
        eigenvalues_[ii][rr] = A[rr][rr];

      if (calculate_eigenvectors) {

        // form groups of equal eigenvalues
        struct Cmp
        {
          bool operator()(const FieldType& a, const FieldType& b) const
          {
            return XT::Common::FloatCmp::lt(a, b);
          }
        };
        std::vector<std::vector<size_t>> eigenvalue_groups;
        std::set<FieldType, Cmp> eigenvalues_done;
        for (size_t jj = 0; jj < rows; ++jj) {
          const auto curr_eigenvalue = eigenvalues_[ii][jj];
          if (!eigenvalues_done.count(curr_eigenvalue)) {
            std::vector<size_t> curr_group;
            curr_group.push_back(jj);
            for (size_t kk = jj + 1; kk < rows; ++kk) {
              if (XT::Common::FloatCmp::eq(curr_eigenvalue, eigenvalues_[ii][kk]))
                curr_group.push_back(kk);
            } // kk
            eigenvalue_groups.push_back(curr_group);
            eigenvalues_done.insert(curr_eigenvalue);
          }
        } // jj

        //         As A Q = Q A_{final} and A_{final} is upper triangular, the first column of Q is always an
        //         eigenvector of A
        //        for (size_t rr = 0; rr < rows; ++rr)
        //          eigenvectors_[ii][rr][0] = (*Q)[rr][0];

        // To get remaining eigenvectors, calculate eigenvectors of A_{final} by solving (A_{final} - \lambda I) x = 0.
        // If x is an eigenvector of A_{final}, Qx is an eigenvector of A.
        for (const auto& group : eigenvalue_groups) {
          size_t value = 1;
          for (const auto& index : group) {
            auto& matrix = *R_k;
            matrix = A;
            for (size_t rr = 0; rr < rows; ++rr)
              matrix[rr][rr] -= eigenvalues_[ii][index];

            VectorType x(0.);
            VectorType rhs = x;
            // backsolve
            for (int rr = int(rows - 1); rr >= 0; --rr) {
              if (XT::Common::FloatCmp::eq(matrix[rr][rr], 0.)) {

                // check if there is a non-zero entry to the right, if so, we can calculate it now
                for (int cc = rr + 1; cc < int(cols); ++cc) {
                  if (XT::Common::FloatCmp::ne(matrix[rr][cc], 0.)) {
                    x[cc] = rhs[rr] / matrix[rr][cc];
                    break;
                  }
                }

                // find first row with nonzero entry and add to all rows above to set entry of this variable to zero
                bool all_rows_zero = true;
                for (int kk = rr - 1; kk >= 0; --kk) {
                  if (XT::Common::FloatCmp::ne(matrix[kk][rr], 0.)) {
                    all_rows_zero = false;

                    // check if there is a non-zero entry to the right, if so, we can set it randomly now
                    for (int cc = rr + 1; cc < int(cols); ++cc) {
                      if (XT::Common::FloatCmp::ne(matrix[kk][cc], 0.)) {
                        x[cc] = FieldType(value++);
                        rhs[kk] -= x[cc] * matrix[kk][cc];
                        matrix[kk][cc] = 0.;
                        break;
                      }
                    }

                    // add current row to rows above
                    for (int ll = kk - 1; ll >= 0; --ll) {
                      const auto factor = -matrix[ll][rr] / matrix[kk][rr];
                      matrix[ll][rr] = 0.;
                      for (int cc = kk; cc < rr; ++cc)
                        matrix[ll][cc] += matrix[kk][cc] * factor;
                      rhs[ll] += rhs[kk] * factor;
                    } // ll
                    break;
                  } // (if mat(kk,rr) != 0)
                } // kk
                if (all_rows_zero)
                  x[rr] = FieldType(value++);
              } else { // if(mat(rr, rr) == 0)

                // check if there is a non-zero entry to the right, if so, we can set it randomly now
                for (int cc = rr + 1; cc < int(cols); ++cc) {
                  if (XT::Common::FloatCmp::ne(matrix[rr][cc], 0.)) {
                    x[cc] = FieldType(value++);
                    rhs[rr] -= x[cc] * matrix[rr][cc];
                    matrix[rr][cc] = 0.;
                    break;
                  }
                }

                x[rr] = rhs[rr] / matrix[rr][rr];

                // set value of x in rows above
                for (int kk = rr - 1; kk >= 0; --kk) {
                  rhs[kk] -= x[rr] * matrix[kk][rr];
                  matrix[kk][rr] = 0.;
                } // kk
              } // else(mat(rr, rr) == 0)
            } // rr

            VectorType Qx(0);
            Q->mv(x, Qx);

            Qx *= 1. / Qx.two_norm();

            for (size_t rr = 0; rr < rows; ++rr)
              *eigenvectors_[ii][rr][index] = Qx[rr];
          } // index

          // orthonormalize eigenvectors in group
          gram_schmidt(ii, group);
        } // groups of eigenvalues
      } // if (calculate_eigenvectors)


      *(eigenvectors_inverse_[ii]) = *(eigenvectors_[ii]);
      eigenvectors_inverse_[ii]->invert();
    } // ii
  }

  void gram_schmidt(const size_t direction, const std::vector<size_t>& indices)
  {
    if (indices.size() > 1) {
      // copy eigenvectors from the matrix eigenvectors_[direction] to a vector of vectors
      std::vector<VectorType> orthonormal_eigenvectors(indices.size());
      for (size_t ii = 0; ii < indices.size(); ++ii)
        for (size_t rr = 0; rr < rows; ++rr)
          orthonormal_eigenvectors[ii][rr] = (*(eigenvectors_[direction]))[rr][indices[ii]];
      // orthonormalize
      for (size_t ii = 1; ii < indices.size(); ++ii) {
        auto& v_i = orthonormal_eigenvectors[ii];
        for (size_t jj = 0; jj < ii; ++jj) {
          const auto& v_j = orthonormal_eigenvectors[jj];
          const auto vj_vj = v_j.dot(v_j);
          const auto vj_vi = v_j.dot(v_i);
          for (size_t rr = 0; rr < rows; ++rr)
            v_i[rr] -= vj_vi / vj_vj * v_j[rr];
        } // jj
        v_i *= 1. / v_i.two_norm();
      } // ii
      // copy eigenvectors back to eigenvectors matrix
      for (size_t ii = 1; ii < indices.size(); ++ii)
        for (size_t rr = 0; rr < rows; ++rr)
          (*(eigenvectors_[direction]))[rr][indices[ii]] = orthonormal_eigenvectors[ii][rr];
    } // if (indices.size() > 1)
  } // void gram_schmidt(...)

  //! \brief modified sign function returning 1 instead of 0 if the value is 0
  FieldType xi(FieldType val) const
  {
    return XT::Common::FloatCmp::eq(val, 0.) ? 1. : val / std::abs(val);
  }

  FieldType get_norm_x(const MatrixType& A, const size_t col_index, size_t num_rows = rows)
  {
    FieldType norm(0);
    for (size_t rr = col_index; rr < num_rows; ++rr)
      norm += std::pow(A[rr][col_index], 2);
    return std::sqrt(norm);
  }

  // Calculates P * A, where P = (I 0 0; 0 I-beta*u*u^T 0; 0 0 I) and u = v[first_row:past_last_row]
  void multiply_householder_from_left(MatrixType& A,
                                      const FieldType& beta,
                                      const VectorType& v,
                                      const size_t first_row = 0,
                                      const size_t past_last_row = rows) const
  {
    // calculate u^T A first
    VectorType uT_A(0.);
    for (size_t cc = 0; cc < cols; ++cc)
      for (size_t rr = first_row; rr < past_last_row; ++rr)
        uT_A[cc] += v[rr] * A[rr][cc];
    // uT_A now contains u^T A[first_row:past_last_row,:]
    for (size_t rr = first_row; rr < past_last_row; ++rr)
      for (size_t cc = 0; cc < cols; ++cc)
        A[rr][cc] -= beta * v[rr] * uT_A[cc];
  }

  // Calculates A * P.
  // \see multiply_householder_from_left
  void multiply_householder_from_right(MatrixType& A,
                                       const FieldType& beta,
                                       const VectorType& v,
                                       const size_t first_col = 0,
                                       const size_t past_last_col = cols) const
  {
    // calculate A u first
    VectorType Au(0.);
    for (size_t rr = 0; rr < rows; ++rr)
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        Au[rr] += A[rr][cc] * v[cc];
    // Au now contains A[:,first_col:past_last_col] u
    for (size_t rr = 0; rr < rows; ++rr)
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        A[rr][cc] -= beta * Au[rr] * v[cc];
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
    R = A;
    for (size_t rr = 0; rr < num_rows; ++rr) {
      Q[rr] *= 0.;
      Q[rr][rr] = 1.;
    }

    VectorType w;
    FieldType tau;
    for (size_t jj = 0; jj < std::min(num_rows - 1, num_cols); ++jj) {
      const auto norm_x = get_norm_x(R, jj, num_rows);
      if (XT::Common::FloatCmp::gt(norm_x, 0.)) {
        // find entry with greatest absolute value for pivoting
        size_t index = jj;
        FieldType max = std::abs(R[jj][jj]);
        for (size_t kk = jj + 1; kk < num_rows; ++kk) {
          if (XT::Common::FloatCmp::gt(std::abs(R[kk][jj]), max)) {
            max = std::abs(R[kk][jj]);
            index = kk;
          }
        }
        if (index != jj) { // swap rows
          // swapping of rows i,j can be done by Householder I - (e_i - e_j)(e_i-e_j)^T
          VectorType e_diff(0.);
          e_diff[jj] = 1.;
          e_diff[index] = -1.;
          multiply_householder_from_left(R, 1., e_diff, jj, num_rows);
          multiply_householder_from_right(Q, 1., e_diff, jj, num_cols);
        }
        const auto s = -sign(R[jj][jj]);
        const FieldType u1 = R[jj][jj] - s * norm_x;
        w[jj] = 1.;
        for (size_t rr = jj + 1; rr < num_rows; ++rr)
          w[rr] = R[rr][jj] / u1;
        tau = -s * u1 / norm_x;

        // calculate R = Q_k R and Q = Q Q_k
        multiply_householder_from_left(R, tau, w, jj, num_rows);
        multiply_householder_from_right(Q, tau, w, jj, num_cols);
      } // if (norm_x != 0)
    } // jj

    // choose Q such that largest entry of each column is positive
    size_t row_index = 0;
    for (size_t cc = 0; cc < num_cols; ++cc) {
      FieldType max = std::numeric_limits<FieldType>::lowest();
      for (size_t rr = 0; rr < num_rows; ++rr) {
        if (std::abs(Q[rr][cc]) > max) {
          max = std::abs(Q[rr][cc]);
          row_index = rr;
        }
      } // rr
      if (XT::Common::FloatCmp::lt(Q[row_index][cc], 0.)) {
        // scal column of Q if largest entry is negative
        // scal row of R to ensure that still A = QR
        for (size_t rr = 0; rr < num_rows; ++rr) {
          Q[rr][cc] = -Q[rr][cc];
          R[cc][rr] = -R[cc][rr];
        } // rr
      } // if (largest entry negative)
    } // cc
  } // void QR_decomp(...)

  //! \brief Transform A to Hessenberg form by transformation P^T A P
  //! \note Expects P to be the unit matrix initially.
  //! \see https://lp.uni-goettingen.de/get/text/2137
  void hessenberg_transformation(MatrixType& A, MatrixType& P) const
  {
    // make P the unit matrix
    for (size_t rr = 0; rr < rows; ++rr) {
      P[rr] *= 0.;
      P[rr][rr] = 1.;
    }
    static_assert(rows == cols, "Hessenberg transformation needs a square matrix!");
    assert(A.N() == rows && A.M() == cols && "A has wrong dimensions!");
    VectorType u(0.);
    for (size_t jj = 0; jj < rows - 2; ++jj) {
      FieldType gamma = 0;
      for (size_t rr = jj + 1; rr < rows; ++rr)
        gamma += std::pow(std::abs(A[rr][jj]), 2);
      gamma = std::sqrt(gamma);
      FieldType beta = gamma * (gamma + std::abs(A[jj + 1][jj]));
      if (XT::Common::FloatCmp::ne(gamma, 0.) && XT::Common::FloatCmp::ne(beta, 0.)) {
        beta = 1. / beta;
        for (size_t rr = jj + 1; rr < rows; ++rr)
          u[rr] = A[rr][jj];
        u[jj + 1] += xi(A[jj + 1][jj]) * gamma;
        // calculate P A P with P = diag(I_j, (I_{n-j} - beta u u*))
        // calculate P A = A - (beta u) (u^T A) first
        multiply_householder_from_left(A, beta, u, jj + 1, rows);
        // now calculate (PA) P  = PA - (PA u) (beta u^T)
        multiply_householder_from_right(A, beta, u, jj + 1, cols);
        // store transformations
        if (calculate_eigenvectors_)
          multiply_householder_from_right(P, beta, u, jj + 1, cols);
      } // if (gamma != 0)
    } // jj
  } // void hessenberg_transformation(...)

  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
  EigenVectorsType eigenvectors_inverse_;
  bool calculate_eigenvectors_;
}; // class QrHouseholderEigenSolver<...>

#if HAVE_LAPACK

template <class FieldType, size_t dimRange>
class UnitMatrix
{
public:
  typedef FieldMatrix<FieldType, dimRange, dimRange> MatrixType;

  // need to reset unit_matrix every time because Lapack changes it
  FieldType* get()
  {
    unit_matrix_ = MatrixType(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      unit_matrix_[rr][rr] = 1.;
    return &(unit_matrix_[0][0]);
  }

private:
  MatrixType unit_matrix_;
}; // class UnitMatrix;

struct LapackWrapper
{
  static int dggev(char jobvl,
                   char jobvr,
                   int n,
                   double* a,
                   int lda,
                   double* b,
                   int ldb,
                   double* alphar,
                   double* alphai,
                   double* beta,
                   double* vl,
                   int ldvl,
                   double* vr,
                   int ldvr);
  static int dggevx(char balanc,
                    char jobvl,
                    char jobvr,
                    char sense,
                    int n,
                    double* a,
                    int lda,
                    double* b,
                    int ldb,
                    double* alphar,
                    double* alphai,
                    double* beta,
                    double* vl,
                    int ldvl,
                    double* vr,
                    int ldvr,
                    int* ilo,
                    int* ihi,
                    double* lscale,
                    double* rscale,
                    double* abnrm,
                    double* bbnrm,
                    double* rconde,
                    double* rcondv);
};

template <class FieldType, size_t dimRange, size_t dimRangeCols>
class LapackEigenSolver
{
  static const size_t rows = dimRange;
  static const size_t cols = dimRange;

public:
  typedef FieldMatrix<FieldType, rows, rows> MatrixType;
  typedef FieldVector<FieldType, rows> VectorType;
  typedef FieldVector<VectorType, dimRangeCols> EigenValuesType;
  typedef FieldVector<std::shared_ptr<MatrixType>, dimRangeCols> EigenVectorsType;
  typedef FieldVector<MatrixType, dimRangeCols> InputMatricesType;

public:
  LapackEigenSolver(InputMatricesType& matrices_in, bool calculate_eigenvectors = false)
  {
    initialize(matrices_in, calculate_eigenvectors);
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
  void initialize(InputMatricesType& matrices_in, const bool calculate_eigenvectors)
  {
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      eigenvectors_[ii] = std::make_shared<MatrixType>();
      eigenvectors_inverse_[ii] = std::make_shared<MatrixType>();

      int N = int(dimRange);
      std::array<double, dimRange> alpha_real, alpha_imag, beta;

      int info = LapackWrapper::dggev('N',
                                      calculate_eigenvectors ? 'V' : 'N',
                                      N,
                                      &(matrices_in[ii][0][0]),
                                      N,
                                      unit_matrix_.get(),
                                      N,
                                      alpha_real.data(),
                                      alpha_imag.data(),
                                      beta.data(),
                                      (double*)nullptr,
                                      N,
                                      &((*(eigenvectors_[ii]))[0][0]),
                                      N);

      //      int ilo, ihi;
      //      std::array<double, dimRange> lscale, rscale, rconde, rcondv;
      //      double abnrm, bbnrm;

      //      int info = LapackWrapper::dggevx('B',
      //                                       'N',
      //                                       calculate_eigenvectors ? 'V' : 'N',
      //                                       'N',
      //                                       N,
      //                                       &(matrices_in[ii][0][0]),
      //                                       N,
      //                                       unit_matrix_.get(),
      //                                       N,
      //                                       alpha_real.data(),
      //                                       alpha_imag.data(),
      //                                       beta.data(),
      //                                       (double*)nullptr,
      //                                       N,
      //                                       &((*(eigenvectors_[ii]))[0][0]),
      //                                       N,
      //                                       &ilo,
      //                                       &ihi,
      //                                       lscale.data(),
      //                                       rscale.data(),
      //                                       &abnrm,
      //                                       &bbnrm,
      //                                       rconde.data(),
      //                                       rcondv.data());

      if (info != 0)
        DUNE_THROW(Dune::MathError, "Lapack returned error " + XT::Common::to_string(info) + "!");

      for (size_t rr = 0; rr < dimRange; ++rr) {
        assert(XT::Common::FloatCmp::eq(alpha_imag[rr], 0.));
        assert(XT::Common::FloatCmp::ne(beta[rr], 0., 1e-6));
        eigenvalues_[ii][rr] = alpha_real[rr] / beta[rr];
      }

      if (calculate_eigenvectors) {
        *(eigenvectors_inverse_[ii]) = *(eigenvectors_[ii]);
        eigenvectors_inverse_[ii]->invert();
      } // if(calculate_eigenvectors)
    } // ii
  } // void initialize(...)

  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
  EigenVectorsType eigenvectors_inverse_;
  static thread_local UnitMatrix<FieldType, dimRange> unit_matrix_;
}; // class LapackEigenSolver<...>

template <class FieldType, size_t dimRange, size_t dimRangeCols>
thread_local UnitMatrix<FieldType, dimRange> LapackEigenSolver<FieldType, dimRange, dimRangeCols>::unit_matrix_;
#endif // HAVE_LAPACK

#if HAVE_EIGEN

template <class FieldType, size_t dimRange, size_t dimRangeCols, bool self_adjoint = false>
class EigenEigenSolver
{
public:
  typedef typename XT::LA::EigenDenseVector<FieldType> VectorType;
  typedef typename XT::LA::EigenDenseMatrix<FieldType> MatrixType;
  typedef FieldVector<FieldVector<FieldType, dimRange>, dimRangeCols> EigenValuesType;
  typedef FieldVector<std::shared_ptr<FieldMatrix<FieldType, dimRange, dimRange>>, dimRangeCols> EigenVectorsType;
  typedef FieldVector<FieldMatrix<FieldType, dimRange, dimRange>, dimRangeCols> InputMatricesType;

private:
  typedef typename MatrixType::BackendType EigenMatrixBackendType;
  typedef typename std::conditional<self_adjoint,
                                    ::Eigen::SelfAdjointEigenSolver<EigenMatrixBackendType>,
                                    ::Eigen::EigenSolver<EigenMatrixBackendType>>::type EigenSolverType;

public:
  EigenEigenSolver(InputMatricesType& matrices_in, bool calculate_eigenvectors = false)
  {
    initialize(matrices_in, calculate_eigenvectors);
  }

  const EigenValuesType& eigenvalues() const
  {
    return eigenvalues_;
  }

  const EigenVectorsType& eigenvectors() const
  {
    return eigenvectors_;
  }

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
  void initialize(InputMatricesType& matrices_in, const bool calculate_eigenvectors)
  {
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      const auto matrix_ii_eigen =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, dimRange>::convert(matrices_in[ii]);
      EigenSolverType eigen_solver(matrix_ii_eigen.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto& eigenvalues_eigen = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      if (XT::Common::FloatCmp::ne(VectorType(eigenvalues_eigen.imag()), VectorType(dimRange, 0.)))
        DUNE_THROW(Dune::MathError, "Eigen returned imaginary eigenvalues!");
      eigenvalues_[ii] = XT::LA::internal::FieldVectorToLaVector<VectorType, dimRange>::convert_back(
          VectorType(eigenvalues_eigen.real()));

      if (calculate_eigenvectors) {
        const auto& eigenvectors_eigen =
            eigen_solver.eigenvectors(); // <- this should be an Eigen vector of std::complex
        if (XT::Common::FloatCmp::ne(MatrixType(eigenvectors_eigen.imag()), MatrixType(dimRange, dimRange, 0.)))
          DUNE_THROW(Dune::MathError, "Eigen returned imaginary eigenvectors!");

        eigenvectors_[ii] = XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, dimRange>::convert_back(
            MatrixType(eigenvectors_eigen.real()));
        eigenvectors_inverse_[ii] =
            XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, dimRange>::convert_back(
                MatrixType(eigenvectors_eigen.real().inverse()));
      } // if (calculate_eigenvectors)
    } // ii
  } // void initialize(...)

  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
  EigenVectorsType eigenvectors_inverse_;
}; // class EigenEigenSolver<...>

#else // HAVE_EIGEN

template <class LocalFluxFunctionImp, bool selfadjoint = false>
class EigenEigenSolver
{
  static_assert(AlwaysFalse<LocalFluxFunctionImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN

template <class FieldType, size_t dimRange, size_t dimRangeCols>
#if HAVE_LAPACK
using DefaultEigenSolver = LapackEigenSolver<FieldType, dimRange, dimRangeCols>;
#elif HAVE_EIGEN
using DefaultEigenSolver = EigenEigenSolver<FieldType, dimRange, dimRangeCols, true>;
#else
using DefaultEigenSolver = QrHouseholderEigenSolver<FieldType, dimRange, dimRangeCols>;
#endif

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
