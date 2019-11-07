// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH

#include <dune/gdt/operators/interfaces.hh>

#include "reconstruction/linear.hh"

namespace Dune {
namespace GDT {


template <class AdvectionOperatorImp, class ReconstructionOperatorImp>
class AdvectionWithReconstructionOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SourceSpaceType::GridViewType,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SourceSpaceType::GridViewType,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;
  using GridViewType = typename AdvectionOperatorType::SourceSpaceType::GridViewType;
  static const size_t r = AdvectionOperatorType::s_r;
  using VectorType = typename AdvectionOperatorType::VectorType;
  using ReconstructionType = DiscreteFunction<VectorType, GridViewType, r>;

  AdvectionWithReconstructionOperator(const AdvectionOperatorType& advection_operator,
                                      const ReconstructionOperatorType& reconstruction_operator)
    : advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
    , reconstruction_(reconstruction_operator.range_space().mapper().size())
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return reconstruction_operator_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    // do reconstruction
    reconstruction_operator_.apply(source, reconstruction_, param);

    // apply advection operator
    std::fill(range.begin(), range.end(), 0.);
    advection_operator_.apply(reconstruction_, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  mutable VectorType reconstruction_;
}; // class AdvectionWithReconstructionOperator<...>


template <class AdvectionOperatorImp, class ReconstructionOperatorImp>
class AdvectionWithPointwiseReconstructionOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SourceSpaceType::GridViewType,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SourceSpaceType::GridViewType,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;
  using GV = typename AdvectionOperatorType::SourceSpaceType::GridViewType;
  static const size_t r = AdvectionOperatorType::s_r;
  using VectorType = typename AdvectionOperatorType::VectorType;

  AdvectionWithPointwiseReconstructionOperator(const AdvectionOperatorType& advection_operator,
                                               const ReconstructionOperatorType& reconstruction_operator)
    : advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
    , reconstructed_values_(source_space().grid_view().indexSet().size(0))
    , reconstructed_function_(source_space().grid_view(), reconstructed_values_)
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
    // do reconstruction
    reconstruction_operator_.apply(source, reconstructed_function_, param);

    // apply advection operator
    std::fill(range.begin(), range.end(), 0.);
    advection_operator_.apply(reconstructed_function_, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  typename ReconstructionOperatorType::ReconstructedValuesType reconstructed_values_;
  mutable typename ReconstructionOperatorType::ReconstructedFunctionType reconstructed_function_;
}; // class AdvectionWithReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_WITH_RECONSTRUCTION_HH
