// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
#define DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH

#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


template <class OperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesOperator
  : public OperatorInterface<typename OperatorImp::MatrixType, typename OperatorImp::SGV, OperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename OperatorImp::MatrixType, typename OperatorImp::SGV, OperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using OperatorType = OperatorImp;
  using InverseHessianOperatorType = InverseHessianOperatorImp;

  EntropicCoordinatesOperator(const OperatorType& operator_in,
                              const InverseHessianOperatorType& inverse_hessian_operator)
    : operator_(operator_in)
    , inverse_hessian_operator_(inverse_hessian_operator)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return operator_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return operator_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    VectorType u_update = range;
    std::fill(u_update.begin(), u_update.end(), 0.);
    operator_.apply(source, u_update, param);
    inverse_hessian_operator_.apply_inverse_hessian(source, u_update, range, param);
  }

  const OperatorType& operator_;
  const InverseHessianOperatorType& inverse_hessian_operator_;
}; // class EntropicCoordinatesOperator<...>

template <class DensityOperatorImp, class AdvectionOperatorImp, class RhsOperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesCombinedOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SGV,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SGV,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using DensityOperatorType = DensityOperatorImp;
  using AdvectionOperatorType = AdvectionOperatorImp;
  using RhsOperatorType = RhsOperatorImp;
  using InverseHessianOperatorType = InverseHessianOperatorImp;

  EntropicCoordinatesCombinedOperator(const DensityOperatorType& density_op,
                                      const AdvectionOperatorType& advection_op,
                                      const RhsOperatorType& rhs_op,
                                      const InverseHessianOperatorType& inverse_hessian_operator)
    : density_op_(density_op)
    , advection_op_(advection_op)
    , rhs_op_(rhs_op)
    , inverse_hessian_operator_(inverse_hessian_operator)
    , reg_indicators_(advection_op_.source_space().grid_view().size(0), false)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return advection_op_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_op_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    density_op_.apply(source, range, param);
    VectorType u_update = range;
    VectorType rhs_update = range;
    std::fill(u_update.begin(), u_update.end(), 0.);
    advection_op_.apply(source, u_update, param);
    u_update *= -1.;
    rhs_op_.apply(source, rhs_update, param);
    u_update += rhs_update;
    std::fill(reg_indicators_.begin(), reg_indicators_.end(), false);
    inverse_hessian_operator_.apply_inverse_hessian(source, u_update, reg_indicators_, range, param);
  }

  const std::vector<bool> reg_indicators() const
  {
    return reg_indicators_;
  }

  const DensityOperatorType& density_op_;
  const AdvectionOperatorType& advection_op_;
  const RhsOperatorType& rhs_op_;
  const InverseHessianOperatorType& inverse_hessian_operator_;
  mutable std::vector<bool> reg_indicators_;
}; // class EntropicCoordinatesOperator<...>


template <class AdvectionOperatorImp, class EntropySolverImp>
class EntropyBasedMomentFvOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SGV,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SGV,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AdvectionOperatorType = AdvectionOperatorImp;
  using EntropySolverType = EntropySolverImp;

  EntropyBasedMomentFvOperator(const AdvectionOperatorType& advection_operator, const EntropySolverType& entropy_solver)
    : advection_operator_(advection_operator)
    , entropy_solver_(entropy_solver)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return advection_operator_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    // solve optimization problems and regularize if necessary
    VectorType regularized = range;
    entropy_solver_.apply(source, regularized, param);

    std::fill(range.begin(), range.end(), 0.);
    advection_operator_.apply(regularized, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const EntropySolverType& entropy_solver_;
}; // class EntropyBasedMomentFvOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
