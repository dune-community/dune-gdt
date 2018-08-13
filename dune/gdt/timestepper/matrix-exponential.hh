// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
#define DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH

#include "config.h"

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>

#include <dune/xt/la/container/common.hh>

#include <dune/gdt/assembler/system.hh>

#include "interface.hh"

#include <dune/xt/data/matrix_exponential/matrix_exponential_extension.hpp>
#include <dune/xt/data/matrix_exponential/matrix_exponential.hpp>

namespace Dune {
namespace GDT {

template <class DiscreteFunctionType, class RhsEvaluationType, class MatrixType>
class MatrixExponentialFunctor
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType> BaseType;
  typedef typename RhsEvaluationType::RangeFieldType FieldType;
  typedef typename RhsEvaluationType::RangeType RangeType;
  static const size_t dimRange = RhsEvaluationType::dimRange;

public:
  using typename BaseType::EntityType;

  MatrixExponentialFunctor(DiscreteFunctionType& solution,
                           const double t,
                           const double dt,
                           const RhsEvaluationType& rhs_evaluation,
                           const std::vector<MatrixType>& matrix_exponentials,
                           const std::vector<MatrixType>& matrix_exponential_integrals)

    : solution_(solution)
    , t_(t)
    , dt_(dt)
    , rhs_evaluation_(rhs_evaluation)
    , matrix_exponentials_(matrix_exponentials)
    , int_exp_mAdt_b_(matrix_exponential_integrals.size())
  {
    for (size_t ii = 0; ii < matrix_exponential_integrals.size(); ++ii) {
      const auto& b = rhs_evaluation_.values()[ii]->b();
      matrix_exponential_integrals[ii].mv(b, int_exp_mAdt_b_[ii]);
    }
  }

  // Solves d_t u(t) = A u(t) + b locally on each entity
  // Multiplying by exp(-At) we get
  // (see https://en.wikipedia.org/wiki/Matrix_exponential#Linear_differential_equations)
  // d_t (exp(-At)u(t)) = exp(-At) b
  // By integrating over (0, dt) wrt. t, we get
  // u(dt) = exp(Adt)(u(0) + (\int_0^{dt} exp(-At)) b)
  virtual void apply_local(const EntityType& entity)
  {
    auto solution_local = solution_.local_discrete_function(entity);

    // get u
    const auto center = entity.geometry().local(entity.geometry().center());
    const auto u0 = solution_local->evaluate(center);

    const size_t subdomain = rhs_evaluation_.subdomain(entity);
    const auto& exp_Adt = matrix_exponentials_[subdomain];

    // calculate solution u = exp(A dt) (u0 + int_exp_mAdt b)
    FieldVector<FieldType, dimRange> ret;
    const auto u = int_exp_mAdt_b_[subdomain] + u0;
    exp_Adt.mv(u, ret);

    // write to return vector
    auto& local_vector = solution_local->vector();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_vector.set(ii, ret[ii]);
  }

private:
  DiscreteFunctionType& solution_;
  const double t_;
  const double dt_;
  const RhsEvaluationType& rhs_evaluation_;
  const std::vector<MatrixType>& matrix_exponentials_;
  std::vector<RangeType> int_exp_mAdt_b_;
};


/** \brief Time stepper solving linear equation d_t u = Au + b by matrix exponential
 */
template <class OperatorImp, class DiscreteFunctionImp>
class MatrixExponentialTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  typedef MatrixExponentialTimeStepper ThisType;
  typedef TimeStepperInterface<DiscreteFunctionImp> BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;
  using typename BaseType::DataHandleType;

  typedef OperatorImp OperatorType;
  typedef typename OperatorType::RhsEvaluationType EvaluationType;
  static const size_t dimDomain = DiscreteFunctionType::dimDomain;
  static const size_t dimRange = DiscreteFunctionType::dimRange;
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> FieldMatrixType;
  typedef XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange> DenseMatrixType;
  typedef XT::LA::CommonSparseOrDenseMatrixCsr<RangeFieldType> SparseMatrixType;

  typedef typename XT::Functions::AffineFluxFunction<typename DiscreteFunctionType::EntityType,
                                                     DomainFieldType,
                                                     dimDomain,
                                                     DiscreteFunctionType,
                                                     RangeFieldType,
                                                     dimRange,
                                                     1>
      RhsAffineFunctionType;
  typedef typename XT::Functions::CheckerboardFunction<typename DiscreteFunctionType::EntityType,
                                                       DomainFieldType,
                                                       dimDomain,
                                                       RangeFieldType,
                                                       dimRange,
                                                       1,
                                                       RhsAffineFunctionType>
      AffineCheckerboardType;


  using BaseType::current_solution;
  using BaseType::current_time;

  MatrixExponentialTimeStepper(const OperatorType& op,
                               DiscreteFunctionType& initial_values,
                               const RangeFieldType t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , op_(op)
    , evaluation_(static_cast<const AffineCheckerboardType&>(op_.evaluation()))
    , num_subdomains_(evaluation_.subdomains())
    , matrix_exponentials_(num_subdomains_, SparseMatrixType(dimRange, dimRange))
    , matrix_exponential_integrals_(num_subdomains_, SparseMatrixType(dimRange, dimRange))
    , last_dt_(0)
  {
  }

  virtual RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    calculate_matrix_exponentials(actual_dt);
    MatrixExponentialFunctor<DiscreteFunctionType, AffineCheckerboardType, SparseMatrixType> functor(
        u_n, t, actual_dt, evaluation_, matrix_exponentials_, matrix_exponential_integrals_);
    SystemAssembler<typename DiscreteFunctionType::SpaceType> assembler(u_n.space());
    assembler.append(functor);
    assembler.assemble(true);

    // augment time
    t += actual_dt;
    last_dt_ = actual_dt;

    return dt;
  } // ... step(...)

private:
  void calculate_matrix_exponentials(const RangeFieldType& actual_dt)
  {
    if (XT::Common::FloatCmp::ne(actual_dt, last_dt_)) {
      size_t num_threads = std::min(XT::Common::threadManager().max_threads(), num_subdomains_);
      std::vector<std::set<size_t>> decomposition(num_threads);
      for (size_t ii = 0; ii < num_subdomains_; ++ii)
        decomposition[ii % num_threads].insert(ii);

      std::vector<std::thread> threads(num_threads);
      // Launch a group of threads
      for (size_t ii = 0; ii < num_threads; ++ii)
        threads[ii] = std::thread(&ThisType::calculate_in_thread, this, actual_dt, decomposition[ii]);
      // Join the threads with the main thread
      for (size_t ii = 0; ii < num_threads; ++ii)
        threads[ii].join();
    }
  }

  void calculate_in_thread(const RangeFieldType& actual_dt, const std::set<size_t>& indices)
  {
    for (const auto& index : indices)
      get_matrix_exponential(index, actual_dt);
  }

  void get_matrix_exponential(size_t index, const RangeFieldType& dt)
  {
    const auto& affine_function = *(evaluation_.values()[index]);
    assert(affine_function.A().size() == 1 && "Not implemented for dimRangeCols > 1!");
    auto A = affine_function.A()[0].operator std::unique_ptr<FieldMatrixType>();

    // calculate matrix exponential exp(A*dt)
    auto Adt = XT::Common::make_unique<DenseMatrixType>(*A);
    *Adt *= dt;
    // get pointer to the underlying array of the FieldMatrix
    RangeFieldType* Adt_array = &((*Adt)[0][0]);
    const double* exp_Adt_array = r8mat_expm1(dimRange, Adt_array);
    auto& exp_Adt = *Adt; // we do not need Adt anymore, use it as storage for exp_Adt
    std::copy_n(exp_Adt_array, dimRange * dimRange, &(exp_Adt[0][0]));
    delete[] exp_Adt_array;
    matrix_exponentials_[index] = SparseMatrixType(exp_Adt, true);

    // calculate \int_0^{dt} exp(-At) dt
    // see https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
    // if A is invertible, the integral is -A^{-1}(exp(-Adt)-I)
    // in general, it is the power series dt*(I + (-Adt)/(2!) + (-Adt)^2/(3!) + ... + (-Adt)^{n-1}/(n!) + ...)
    auto A_inverse = XT::Common::make_unique<DenseMatrixType>(*A);
    auto& int_exp_mAdt = *Adt; // use Adt's storage again
    try {
      // For 1x1, 2x2 and 3x3 matrices, dune-common only checks whether the matrix is actually invertible if
      // DUNE_FMatrix_WITH_CHECKING is set. We always want to check to
      // be able to catch the error, so we have to copy the checks over from dune/common/fmatrix.hh.
      if (dimRange == 1) {
        if (fvmeta::absreal((*A)[0][0]) < FMatrixPrecision<>::absolute_limit())
          DUNE_THROW(FMatrixError, "matrix is singular");
      } else if (dimRange == 2) {
        RangeFieldType det = (*A)[0][0] * (*A)[1][1] - (*A)[0][1] * (*A)[1][0];
        if (fvmeta::absreal(det) < FMatrixPrecision<>::absolute_limit())
          DUNE_THROW(FMatrixError, "matrix is singular");
      } else if (dimRange == 3) {
        RangeFieldType t4 = (*A)[0][0] * (*A)[1][1];
        RangeFieldType t6 = (*A)[0][0] * (*A)[1][2];
        RangeFieldType t8 = (*A)[0][1] * (*A)[1][0];
        RangeFieldType t10 = (*A)[0][2] * (*A)[1][0];
        RangeFieldType t12 = (*A)[0][1] * (*A)[2][0];
        RangeFieldType t14 = (*A)[0][2] * (*A)[2][0];
        RangeFieldType det = (t4 * (*A)[2][2] - t6 * (*A)[2][1] - t8 * (*A)[2][2] + t10 * (*A)[2][1] + t12 * (*A)[1][2]
                              - t14 * (*A)[1][1]);
        if (fvmeta::absreal(det) < FMatrixPrecision<>::absolute_limit())
          DUNE_THROW(FMatrixError, "matrix is singular");
      }
      A_inverse->invert();
      *A_inverse *= -1.;

      // calculate matrix exponential exp(-A*dt)
      auto* mAdt = A.get();
      *mAdt *= -dt;
      // get pointer to the underlying array of the FieldMatrix
      double* mAdt_array = &((*mAdt)[0][0]);
      const double* exp_mAdt_array = r8mat_expm1(dimRange, mAdt_array);
      std::copy_n(exp_mAdt_array, dimRange * dimRange, &(int_exp_mAdt[0][0]));
      delete[] exp_mAdt_array;
      for (size_t ii = 0; ii < dimRange; ++ii)
        int_exp_mAdt[ii][ii] -= 1.;
      int_exp_mAdt.leftmultiply(*A_inverse);
    } catch (Dune::FMatrixError&) { // A not invertible
      auto* minus_A = A.get();
      *minus_A *= -1.;
      const double* int_exp_mAdt_array = r8mat_expm_integral(dimRange, &((*minus_A)[0][0]), dt);
      std::copy_n(int_exp_mAdt_array, dimRange * dimRange, &(int_exp_mAdt[0][0]));
      delete[] int_exp_mAdt_array;
    } // catch(...)
    matrix_exponential_integrals_[index] = SparseMatrixType(int_exp_mAdt, true);
  } // void get_matrix_exponential(...)

  const OperatorType& op_;
  const AffineCheckerboardType& evaluation_;
  size_t num_subdomains_;
  std::vector<SparseMatrixType> matrix_exponentials_;
  std::vector<SparseMatrixType> matrix_exponential_integrals_;
  RangeFieldType last_dt_;
};

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
