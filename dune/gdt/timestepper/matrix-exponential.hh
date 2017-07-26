// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
#define DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/checkerboard.hh>

#include <dune/xt/la/container/common.hh>

#include <dune/gdt/assembler/system.hh>

#include "interface.hh"

#include "matrixexponential/matrix_exponential.hpp"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType, class RhsEvaluationType>
class MatrixExponentialFunctor
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType> BaseType;
  typedef typename RhsEvaluationType::RangeFieldType FieldType;
  static const size_t dimRange = RhsEvaluationType::dimRange;

public:
  using typename BaseType::EntityType;
  typedef FieldMatrix<FieldType, dimRange, dimRange> MatrixType;

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
    , matrix_exponential_integrals_(matrix_exponential_integrals)
  {
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
    const auto& b = rhs_evaluation_.values()[subdomain]->b();
    const auto& exp_Adt = matrix_exponentials_[subdomain];
    const auto& int_exp_mAdt = matrix_exponential_integrals_[subdomain];

    // calculate solution u = exp(A dt) (u0 + int_exp_mAdt b)
    FieldVector<FieldType, dimRange> u(0), ret;
    int_exp_mAdt.mv(b, u);
    u += u0;
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
  const std::vector<MatrixType>& matrix_exponential_integrals_;
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
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;

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
                               const DiscreteFunctionType& initial_values,
                               const RangeFieldType t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , op_(op)
    , evaluation_(static_cast<const AffineCheckerboardType&>(op_.evaluation()))
    , num_subdomains_(evaluation_.subdomains())
    , matrix_exponentials_(num_subdomains_)
    , matrix_exponential_integrals_(num_subdomains_)
    , last_dt_(0)
  {
  }

  virtual RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    calculate_matrix_exponentials(actual_dt);
    MatrixExponentialFunctor<DiscreteFunctionType, AffineCheckerboardType> functor(
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
    auto A = affine_function.A()[0].operator std::unique_ptr<MatrixType>();

    // calculate matrix exponential exp(A*dt)
    auto Adt = XT::Common::make_unique<MatrixType>(*A);
    *Adt *= dt;
    // get pointer to the underlying array of the FieldMatrix
    RangeFieldType* Adt_array = &((*Adt)[0][0]);
    const double* exp_Adt_array = r8mat_expm1(dimRange, Adt_array);
    auto& exp_Adt = matrix_exponentials_[index];
    std::copy_n(exp_Adt_array, dimRange * dimRange, &(exp_Adt[0][0]));
    delete[] exp_Adt_array;

    // calculate \int_0^{dt} exp(-At) dt
    // see https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
    // if A is invertible, the integral is -A^{-1}(exp(-Adt)-I)
    // in general, it is the power series dt*(I + (-Adt)/(2!) + (-Adt)^2/(3!) + ... + (-Adt)^{n-1}/(n!) + ...)
    auto& int_exp_mAdt = matrix_exponential_integrals_[index];
    auto A_inverse = XT::Common::make_unique<MatrixType>(*A);
    try {
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
  } // void get_matrix_exponential(...)

  const OperatorType& op_;
  const AffineCheckerboardType& evaluation_;
  size_t num_subdomains_;
  std::vector<MatrixType> matrix_exponentials_;
  std::vector<MatrixType> matrix_exponential_integrals_;
  RangeFieldType last_dt_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
