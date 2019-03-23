// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH

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

  virtual bool linear() const override final
  {
    return false;
  }

  virtual const SourceSpaceType& source_space() const override final
  {
    return reconstruction_operator_.source_space();
  }

  virtual const RangeSpaceType& range_space() const override final
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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH
