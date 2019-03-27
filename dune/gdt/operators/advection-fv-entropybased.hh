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

  virtual bool linear() const override final
  {
    return false;
  }

  virtual const SourceSpaceType& source_space() const override final
  {
    return advection_operator_.source_space();
  }

  virtual const RangeSpaceType& range_space() const override final
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
