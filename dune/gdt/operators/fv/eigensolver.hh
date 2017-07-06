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

extern "C" void dggev_(const char* JOBVL,
                       const char* JOBVR,
                       const int* N,
                       const double* A,
                       const int* LDA,
                       const double* B,
                       const int* LDB,
                       double* ALPHAR,
                       double* ALPHAI,
                       double* BETA,
                       double* VL,
                       const int* LDVL,
                       double* VR,
                       const int* LDVR,
                       double* WORK,
                       const int* LWORK,
                       int* INFO);


namespace Dune {
namespace GDT {

template <class LocalFluxFunctionImp>
class LapackEigenSolver
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
  typedef FieldVector<FieldVector<RangeFieldType, dimRange>, dimRangeCols> EigenValuesType;
  typedef FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimRangeCols> EigenVectorsType;

private:
  typedef typename MatrixType::BackendType EigenMatrixBackendType;

public:
  LapackEigenSolver(const LocalFluxFunctionImp& local_flux_function,
                    const DomainType& x_local,
                    const StateRangeType& u,
                    const XT::Common::Parameter& param)
  {
    const auto partial_u = XT::Functions::JacobianRangeTypeConverter<stateDimRange, dimRange, dimRangeCols>::convert(
        local_flux_function.partial_u(x_local, u, param));
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {

      int N = int(dimRange);
      double alpha_real[dimRange], alpha_imag[dimRange], beta[dimRange];
      double workdummy;
      int lwork = -1; // Request optimum work size.
      int info = 0;

      FieldMatrix<RangeFieldType, dimRange, dimRange> I(0.);
      for (size_t rr = 0; rr < dimRange; ++rr)
        I[rr][rr] = 1.;

      // lapack uses column-major representation, so transpose first
      auto A = partial_u[ii];
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < dimRange; ++cc)
          A[cc][rr] = partial_u[ii][rr][cc];

      ::dggev_("N",
               "V",
               &N,
               &(A[0][0]),
               &N,
               &(I[0][0]),
               &N,
               alpha_real,
               alpha_imag,
               beta,
               nullptr,
               &N,
               &(eigenvectors_[ii][0][0]),
               &N,
               &workdummy,
               &lwork,
               &info);


      lwork = int(workdummy + 0.5);
      std::vector<double> work(lwork);

      ::dggev_("N",
               "V",
               &N,
               &(A[0][0]),
               &N,
               &(I[0][0]),
               &N,
               alpha_real,
               alpha_imag,
               beta,
               nullptr,
               &N,
               &(eigenvectors_[ii][0][0]),
               &N,
               work.data(),
               &lwork,
               &info);

      // transpose eigenvectors, as lapack uses column-major
      auto vec_copy = eigenvectors_[ii];
      for (size_t rr = 0; rr < dimRange; ++rr)
        for (size_t cc = 0; cc < dimRange; ++cc)
          eigenvectors_[ii][cc][rr] = vec_copy[rr][cc];

      for (size_t rr = 0; rr < dimRange; ++rr) {
        assert(XT::Common::FloatCmp::eq(alpha_imag[rr], 0.));
        assert(XT::Common::FloatCmp::eq(beta[rr], 1.));
        eigenvalues_[ii][rr] = alpha_real[rr];
      }

      eigenvectors_inverse_[ii] = eigenvectors_[ii];
      eigenvectors_inverse_[ii].invert();
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

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
  EigenValuesType eigenvalues_;
  EigenVectorsType eigenvectors_;
  EigenVectorsType eigenvectors_inverse_;
}; // class LapackEigenSolver<...>

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
  typedef FieldVector<FieldVector<RangeFieldType, dimRange>, dimRangeCols> EigenValuesType;
  typedef FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimRangeCols> EigenVectorsType;

public:
  QrHouseholderEigenSolver(const LocalFluxFunctionImp& local_flux_function,
                           const DomainType& x_local,
                           const StateRangeType& u,
                           const XT::Common::Parameter& param,
                           const bool calculate_eigenvectors = true)
    : tmp_vec_(u.size())
    , calculate_eigenvectors_(calculate_eigenvectors)
  {
    const auto partial_u = XT::Functions::JacobianRangeTypeConverter<stateDimRange, dimRange, dimRangeCols>::convert(
        local_flux_function.partial_u(x_local, u, param));

    static constexpr size_t max_counts = 10000;
    MatrixType A, R_k, Q_k(rows, rows), Q(rows, rows), Q_copy(rows, rows);
    const RangeFieldType tol = 1e-15;
    for (size_t ii = 0; ii < dimRangeCols; ++ii) {
      for (size_t rr = 0; rr < rows; ++rr)
        Q.unit_row(rr);
      A = XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, stateDimRange>::convert(partial_u[ii]);
      hessenberg_transformation(A, Q);
      std::cout << "A: " << XT::Common::to_string(A, 20) << std::endl;
      std::cout << "Q: " << XT::Common::to_string(Q, 20) << std::endl;
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

          std::cout << "A: " << XT::Common::to_string(A, 20) << std::endl;
          std::cout << "Q_k: " << XT::Common::to_string(Q_k, 20) << std::endl;
          std::cout << "R_k: " << XT::Common::to_string(R_k, 20) << std::endl;

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

          auto A_copy = A;
          // update upper right part
          for (size_t rr = 0; rr < num_rows; ++rr) {
            for (size_t cc = num_cols; cc < cols; ++cc) {
              A.set_entry(rr, cc, 0.);
              for (size_t ll = 0; ll < num_rows; ++ll)
                A.add_to_entry(rr, cc, Q_k.get_entry(ll, rr) * A_copy.get_entry(ll, cc));
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
        eigenvalues_[ii][rr] = A.get_entry(rr, rr);

      if (calculate_eigenvectors_) {

        // form groups of equal eigenvalues
        struct Cmp
        {
          bool operator()(const RangeFieldType& a, const RangeFieldType& b) const
          {
            return XT::Common::FloatCmp::lt(a, b);
          }
        };
        std::vector<std::vector<size_t>> eigenvalue_groups;
        std::set<RangeFieldType, Cmp> eigenvalues_done;
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

        // As A Q = Q A_{final} and A_{final} is upper triangular, the first column of Q is always an eigenvector of A
        for (size_t rr = 0; rr < rows; ++rr)
          eigenvectors_[ii][rr][0] = Q.get_entry(rr, 0);

        // To get remaining eigenvectors, calculate eigenvectors of A_{final} by solving (A_{final} - \lambda I) x = 0.
        // If x is an eigenvector of A_{final}, Qx is an eigenvector of A.
        for (const auto& group : eigenvalue_groups) {
          size_t value = 1;
          for (const auto& index : group) {
            auto matrix = A.copy();
            for (size_t rr = 0; rr < rows; ++rr)
              matrix.add_to_entry(rr, rr, -eigenvalues_[ii][index]);

            VectorType x(rows, 0.);
            VectorType rhs = x;
            // backsolve
            for (int rr = int(rows - 1); rr >= 0; --rr) {
              if (XT::Common::FloatCmp::eq(matrix.get_entry(rr, rr), 0.)) {

                // check if there is a non-zero entry to the right, if so, we can calculate it now
                for (int cc = rr + 1; cc < int(cols); ++cc) {
                  if (XT::Common::FloatCmp::ne(matrix.get_entry(rr, cc), 0.)) {
                    if (std::abs(rhs.get_entry(rr) / matrix.get_entry(rr, cc) > 1000))
                      std::cout << rhs.get_entry(rr) / matrix.get_entry(rr, cc) << " " << matrix.get_entry(rr, cc)
                                << " " << rhs.get_entry(rr) << std::endl;
                    x.set_entry(cc, rhs.get_entry(rr) / matrix.get_entry(rr, cc));
                    break;
                  }
                }

                // find first row with nonzero entry and add to all rows above to set entry of this variable to zero
                bool all_rows_zero = true;
                for (int kk = rr - 1; kk >= 0; --kk) {
                  if (XT::Common::FloatCmp::ne(matrix.get_entry(kk, rr), 0.)) {
                    all_rows_zero = false;

                    // check if there is a non-zero entry to the right, if so, we can set it randomly now
                    for (int cc = rr + 1; cc < int(cols); ++cc) {
                      if (XT::Common::FloatCmp::ne(matrix.get_entry(kk, cc), 0.)) {
                        x.set_entry(cc, RangeFieldType(value++));
                        rhs.add_to_entry(kk, -x.get_entry(cc) * matrix.get_entry(kk, cc));
                        matrix.set_entry(kk, cc, 0.);
                        break;
                      }
                    }

                    // add current row to rows above
                    for (int ll = kk - 1; ll >= 0; --ll) {
                      const auto factor = -matrix.get_entry(ll, rr) / matrix.get_entry(kk, rr);
                      if (std::abs(factor) > 1000)
                        std::cout << factor << " " << kk << " " << ll << " " << rr << " " << matrix.get_entry(ll, rr)
                                  << " " << matrix.get_entry(kk, rr) << std::endl;
                      matrix.set_entry(ll, rr, 0.);
                      for (int cc = kk; cc < rr; ++cc)
                        matrix.add_to_entry(ll, cc, matrix.get_entry(kk, cc) * factor);
                      rhs.add_to_entry(ll, rhs.get_entry(kk) * factor);
                    } // ll
                    break;
                  } // (if mat(kk,rr) != 0)
                } // kk
                if (all_rows_zero)
                  x.set_entry(rr, RangeFieldType(value++));
              } else { // if(mat(rr, rr) == 0)

                // check if there is a non-zero entry to the right, if so, we can set it randomly now
                for (int cc = rr + 1; cc < int(cols); ++cc) {
                  if (XT::Common::FloatCmp::ne(matrix.get_entry(rr, cc), 0.)) {
                    x.set_entry(cc, RangeFieldType(value++));
                    rhs.add_to_entry(rr, -x.get_entry(cc) * matrix.get_entry(rr, cc));
                    matrix.set_entry(rr, cc, 0.);
                    break;
                  }
                }

                x.set_entry(rr, rhs.get_entry(rr) / matrix.get_entry(rr, rr));

                // set value of x in rows above
                for (int kk = rr - 1; kk >= 0; --kk) {
                  rhs.add_to_entry(kk, -x.get_entry(rr) * matrix.get_entry(kk, rr));
                  matrix.set_entry(kk, rr, 0.);
                } // kk
              } // else(mat(rr, rr) == 0)
            } // rr

            std::cout << "A: " << XT::Common::to_string(A, 20) << std::endl;
            std::cout << "Q: " << XT::Common::to_string(Q, 20) << std::endl;
            //            std::cout << "matrix: " << XT::Common::to_string(matrix, 20) << std::endl;
            std::cout << "eigvec: " << XT::Common::to_string(x, 20) << std::endl;
            std::cout << "eigval: " << XT::Common::to_string(eigenvalues_[ii][index], 20) << std::endl;
            Q.mv(x, *tmp_vec_);
            std::cout << "eigvec after: " << XT::Common::to_string(*tmp_vec_, 20) << std::endl;

            std::cout << "l2_norm() " << tmp_vec_->l2_norm() << std::endl;
            tmp_vec_->scal(1. / tmp_vec_->l2_norm());
            std::cout << "eigvec after scaling: " << XT::Common::to_string(*tmp_vec_, 20) << std::endl;

            for (size_t rr = 0; rr < rows; ++rr)
              eigenvectors_[ii][rr][index] = tmp_vec_->get_entry(rr);
          } // index

          // orthonormalize eigenvectors in group
          gram_schmidt(ii, group);
        } // groups of eigenvalues
      } // if (calculate_eigenvectors_)

      //      XT::LA::Solver<MatrixType> solver(eigenvectors_[ii]);
      //      for (size_t cc = 0; cc < cols; ++cc) {
      //         as A A^{-1} = I, solve A a_inv_j = e_j where a_inv_j is the j-th column of A^{-1}
      //        VectorType rhs(rows, 0);
      //        rhs.set_entry(cc, 1.);
      //        solver.apply(rhs, *tmp_vec_);
      //        for (size_t rr = 0; rr < rows; ++rr)
      //          eigenvectors_inverse_[ii].set_entry(rr, cc, tmp_vec_->get_entry(rr));
      //      } // cc
      std::cout << "A: " << XT::Common::to_string(partial_u[ii], 20) << std::endl;
      std::cout << "eigvals: " << XT::Common::to_string(eigenvalues_[ii], 20) << std::endl;
      std::cout << "eigvecs: " << XT::Common::to_string(eigenvectors_[ii], 20) << std::endl;
      eigenvectors_inverse_[ii] = eigenvectors_[ii];
      eigenvectors_inverse_[ii].invert();
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

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
  void gram_schmidt(const size_t direction, const std::vector<size_t>& indices)
  {
    if (indices.size() > 1) {
      // copy eigenvectors from the matrix eigenvectors_[direction] to a vector of vectors
      std::vector<VectorType> orthonormal_eigenvectors(indices.size(), VectorType(rows));
      for (size_t ii = 0; ii < indices.size(); ++ii)
        for (size_t rr = 0; rr < rows; ++rr)
          orthonormal_eigenvectors[ii].set_entry(rr, eigenvectors_[direction][rr][indices[ii]]);
      // orthonormalize
      for (size_t ii = 1; ii < indices.size(); ++ii) {
        auto& v_i = orthonormal_eigenvectors[ii];
        for (size_t jj = 0; jj < ii; ++jj) {
          const auto& v_j = orthonormal_eigenvectors[jj];
          const auto vj_vj = v_j.dot(v_j);
          const auto vj_vi = v_j.dot(v_i);
          for (size_t rr = 0; rr < rows; ++rr)
            v_i.add_to_entry(rr, -vj_vi / vj_vj * v_j.get_entry(rr));
        } // jj
        v_i.scal(1. / v_i.l2_norm());
      } // ii
      // copy eigenvectors back to eigenvectors matrix
      for (size_t ii = 1; ii < indices.size(); ++ii)
        for (size_t rr = 0; rr < rows; ++rr)
          eigenvectors_[direction][rr][indices[ii]] = orthonormal_eigenvectors[ii].get_entry(rr);
    } // if (indices.size() > 1)
  } // void gram_schmidt(...)

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

  // Calculates P * A, where P = (I 0 0; 0 I-beta*u*u^T 0; 0 0 I) and u = v[first_row:past_last_row]
  void multiply_householder_from_left(MatrixType& A,
                                      const RangeFieldType& beta,
                                      const VectorType& v,
                                      const size_t first_row = 0,
                                      const size_t past_last_row = rows) const
  {
    // calculate u^T A first
    for (size_t cc = 0; cc < cols; ++cc) {
      tmp_vec_->set_entry(cc, 0.);
      for (size_t rr = first_row; rr < past_last_row; ++rr)
        tmp_vec_->add_to_entry(cc, v.get_entry(rr) * A.get_entry(rr, cc));
    } // tmp_vec now contains u^T A[first_row:past_last_row,:]
    for (size_t rr = first_row; rr < past_last_row; ++rr)
      for (size_t cc = 0; cc < cols; ++cc)
        A.add_to_entry(rr, cc, -beta * v.get_entry(rr) * tmp_vec_->get_entry(cc));
  }

  // Calculates A * P.
  // \see multiply_householder_from_left
  void multiply_householder_from_right(MatrixType& A,
                                       const RangeFieldType& beta,
                                       const VectorType& v,
                                       const size_t first_col = 0,
                                       const size_t past_last_col = cols) const
  {
    // calculate A u first
    for (size_t rr = 0; rr < rows; ++rr) {
      tmp_vec_->set_entry(rr, 0.);
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        tmp_vec_->add_to_entry(rr, A.get_entry(rr, cc) * v.get_entry(cc));
    } // tmp_vec now contains A[:,first_col:past_last_col] u
    for (size_t rr = 0; rr < rows; ++rr)
      for (size_t cc = first_col; cc < past_last_col; ++cc)
        A.add_to_entry(rr, cc, -beta * tmp_vec_->get_entry(rr) * v.get_entry(cc));
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

    VectorType w(num_rows);
    RangeFieldType tau;
    for (size_t jj = 0; jj < std::min(num_rows - 1, num_cols); ++jj) {
      const auto norm_x = get_norm_x(R, jj, num_rows);
      if (XT::Common::FloatCmp::gt(norm_x, 0.)) {
#if 0
        // find entry with greatest absolute value
        size_t index = jj;
        RangeFieldType max = std::abs(R.get_entry(jj, jj));
        for (size_t kk = jj + 1; kk < num_rows; ++kk) {
          if (XT::Common::FloatCmp::gt(std::abs(R.get_entry(kk, jj)), max)) {
            max = std::abs(R.get_entry(kk, jj));
            index = kk;
          }
        }
        if (index != jj) { // swap rows
          // swapping of rows i,j can be done by Householder I - (e_i - e_j)(e_i-e_j)^T
          VectorType e_diff(num_rows, 0.);
          e_diff.set_entry(jj, 1.);
          e_diff.set_entry(index, -1.);
          multiply_householder_from_left(R, 1., e_diff, jj, num_rows);
          multiply_householder_from_right(Q, 1., e_diff, jj, num_cols);
        }
#endif
        const auto s = -sign(R.get_entry(jj, jj));
        const RangeFieldType u1 = R.get_entry(jj, jj) - s * norm_x;
        w.set_entry(jj, 1.);
        for (size_t rr = jj + 1; rr < num_rows; ++rr)
          w.set_entry(rr, R.get_entry(rr, jj) / u1);
        tau = -s * u1 / norm_x;

        // calculate R = Q_k R and Q = Q Q_k
        multiply_householder_from_left(R, tau, w, jj, num_rows);
        multiply_householder_from_right(Q, tau, w, jj, num_cols);
      } // if (norm_x != 0)
    } // jj

    // choose Q such that largest entry of each column is positive
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

  //! \brief Transform A to Hessenberg form by transformation P^T A P
  //! \note Expects P to be the unit matrix initially.
  //! \see https://lp.uni-goettingen.de/get/text/2137
  void hessenberg_transformation(MatrixType& A, MatrixType& P) const
  {
    static_assert(rows == cols, "Hessenberg transformation needs a square matrix!");
    assert(A.rows() == rows && A.cols() == cols && "A has wrong dimensions!");
    VectorType u(rows);
    for (size_t jj = 0; jj < rows - 2; ++jj) {
      RangeFieldType gamma = 0;
      for (size_t rr = jj + 1; rr < rows; ++rr)
        gamma += std::pow(std::abs(A.get_entry(rr, jj)), 2);
      gamma = std::sqrt(gamma);
      RangeFieldType beta = gamma * (gamma + std::abs(A.get_entry(jj + 1, jj)));
      if (XT::Common::FloatCmp::ne(gamma, 0.) && XT::Common::FloatCmp::ne(beta, 0.)) {
        beta = 1. / beta;
        for (size_t rr = jj + 1; rr < rows; ++rr)
          u.set_entry(rr, A.get_entry(rr, jj));
        u.add_to_entry(jj + 1, xi(A.get_entry(jj + 1, jj)) * gamma);
        std::cout << "beta: " << beta << ", gamma: " << gamma << ", Aentry: " << A.get_entry(jj + 1, jj) << std::endl;
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
  mutable XT::Common::PerThreadValue<VectorType> tmp_vec_;
  const bool calculate_eigenvectors_;
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
  typedef FieldVector<FieldVector<RangeFieldType, dimRange>, dimRangeCols> EigenValuesType;
  typedef FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimRangeCols> EigenVectorsType;

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
      std::cout << "A1: " << XT::Common::to_string(partial_u[ii]) << std::endl;
      std::cout << "A2: " << XT::Common::to_string(partial_u_eigen) << std::endl;
      EigenSolverType eigen_solver(partial_u_eigen.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto& eigenvalues_eigen = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto& eigenvectors_eigen = eigen_solver.eigenvectors(); // <- this should be an Eigen vector of std::complex
      assert(size_t(eigenvalues_eigen.size()) == dimRange);
      std::cout << "real: " << XT::Common::to_string(VectorType(eigenvalues_eigen.real())) << std::endl;
      std::cout << "imag: " << XT::Common::to_string(VectorType(eigenvalues_eigen.imag())) << std::endl;
      std::cout << "real vec: " << XT::Common::to_string(MatrixType(eigenvectors_eigen.real())) << std::endl;
      std::cout << "imag vec: " << XT::Common::to_string(MatrixType(eigenvectors_eigen.imag())) << std::endl;
      assert(XT::Common::FloatCmp::eq(VectorType(eigenvalues_eigen.imag()), VectorType(dimRange, 0.)));
      assert(XT::Common::FloatCmp::eq(MatrixType(eigenvectors_eigen.imag()), MatrixType(dimRange, dimRange, 0.)));
      eigenvalues_[ii] = XT::LA::internal::FieldVectorToLaVector<VectorType, dimRange>::convert_back(
          VectorType(eigenvalues_eigen.real()));
      eigenvectors_[ii] =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, stateDimRange>::convert_back(
              MatrixType(eigenvectors_eigen.real()));
      eigenvectors_inverse_[ii] =
          XT::LA::internal::FieldMatrixToLaDenseMatrix<MatrixType, dimRange, stateDimRange>::convert_back(
              MatrixType(eigenvectors_eigen.real().inverse()));
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

  const EigenVectorsType& eigenvectors_inverse() const
  {
    return eigenvectors_inverse_;
  }

private:
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

template <class LocalFluxFunctionImp>
//#if HAVE_EIGEN
// using DefaultEigenSolver = EigenEigenSolver<LocalFluxFunctionImp>;
//#else
// using DefaultEigenSolver = QrHouseholderEigenSolver<LocalFluxFunctionImp>;
//#endif
using DefaultEigenSolver = LapackEigenSolver<LocalFluxFunctionImp>;

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_EIGENSOLVER_HH
