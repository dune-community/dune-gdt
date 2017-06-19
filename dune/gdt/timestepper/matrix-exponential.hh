// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
#define DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH

#include <dune/gdt/assembler/system.hh>

#include "interface.hh"

#include "matrixexponential/matrix_exponential.hpp"

namespace Dune {
namespace GDT {


// TODO: make thread-safe
template <class DiscreteFunctionType, class RhsEvaluationType>
class MatrixExponentialFunctor
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType> BaseType;
  typedef typename RhsEvaluationType::RangeFieldType FieldType;
  static const size_t dimRange = RhsEvaluationType::dimRange;

public:
  using typename BaseType::EntityType;

  MatrixExponentialFunctor(DiscreteFunctionType& solution,
                           const double t,
                           const double dt,
                           const RhsEvaluationType& rhs_evaluation)
    : solution_(solution)
    , t_(t)
    , dt_(dt)
    , rhs_evaluation_(rhs_evaluation)
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

    // get A and b
    auto zero = u0;
    zero *= 0.;
    const auto local_rhs = rhs_evaluation_.local_function(entity);
    FieldMatrix<FieldType, dimRange, dimRange> A = local_rhs->jacobian_wrt_u(center, u0);
    const auto b = local_rhs->evaluate(center, u0);

    // calculate matrix exponential exp(A*dt)
    auto Adt = A;
    Adt *= dt_;
    // get pointer to the underlying array of the FieldMatrix
    double* Adt_array = &(Adt[0][0]);

    const double* exp_Adt_array = r8mat_expm1(dimRange, Adt_array);

    FieldMatrix<FieldType, dimRange, dimRange> exp_Adt;
    std::copy_n(exp_Adt_array, dimRange * dimRange, &(exp_Adt[0][0]));
    delete[] exp_Adt_array;

    // calculate integral of exp(-At) int_exp_mAt
    // see https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
    // if A is invertible, the integral is -A^{-1}(exp(-Adt)-I)
    // in general, it is the power series dt*(I + (-Adt)/(2!) + (-Adt)^2/(3!) + ... + (-Adt)^{n-1}/(n!) + ...)
    FieldMatrix<FieldType, dimRange, dimRange> int_exp_mAdt;
    try {
      auto A_inverse = A;
      A_inverse.invert();
      A_inverse *= -1.;

      // calculate matrix exponential exp(-A*dt)
      auto mAdt = A;
      mAdt *= -dt_;
      // get pointer to the underlying array of the FieldMatrix
      double* mAdt_array = &(mAdt[0][0]);
      const double* exp_mAdt_array = r8mat_expm1(dimRange, mAdt_array);

      std::copy_n(exp_mAdt_array, dimRange * dimRange, &(int_exp_mAdt[0][0]));
      delete[] exp_mAdt_array;

      for (size_t ii = 0; ii < dimRange; ++ii)
        int_exp_mAdt[ii][ii] -= 1.;
      int_exp_mAdt.leftmultiply(A_inverse);
    } catch (Dune::FMatrixError&) {
      auto minus_A = A;
      minus_A *= -1.;
      const double* int_exp_mAdt_array = r8mat_expm_integral(dimRange, &(minus_A[0][0]), dt_);
      std::copy_n(int_exp_mAdt_array, dimRange * dimRange, &(int_exp_mAdt[0][0]));
    }

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
};


/** \brief Time stepper solving linear equation d_t u = Au + b by matrix exponential
 */
template <class OperatorImp, class DiscreteFunctionImp, class TimeFieldImp = double>
class MatrixExponentialTimeStepper : public TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp> BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;
  using typename BaseType::DataHandleType;

  typedef OperatorImp OperatorType;

  using BaseType::current_solution;
  using BaseType::current_time;

  MatrixExponentialTimeStepper(const OperatorType& op,
                               const DiscreteFunctionType& initial_values,
                               const TimeFieldImp t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , op_(op)
  {
  }

  virtual TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt) override final
  {
    const TimeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();

    MatrixExponentialFunctor<DiscreteFunctionType, typename OperatorType::RhsEvaluationType> functor(
        u_n, t, actual_dt, op_.evaluation());
    SystemAssembler<typename DiscreteFunctionType::SpaceType> assembler(u_n.space());
    assembler.append(functor);
    assembler.assemble(true);

    // augment time
    t += actual_dt;

    return dt;
  } // ... step(...)

private:
  const OperatorType& op_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_MATRIXEXPONENTIAL_HH
