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


template <class AdvectionOperatorImp, class ReconstructionOperatorImp, class Traits>
class AdvectionWithReconstructionOperator;


namespace internal {


template <class AdvectionOperatorImp, class ReconstructionOperatorImp>
class AdvectionWithReconstructionOperatorTraits
{
public:
  using AdvectionOperatorType = AdvectionOperatorImp;
  using FieldType = typename AdvectionOperatorType::FieldType;
  using JacobianType = NoJacobian;
  using ReconstructionOperatorType = ReconstructionOperatorImp;
  using derived_type = AdvectionWithReconstructionOperator<AdvectionOperatorImp,
                                                           ReconstructionOperatorImp,
                                                           AdvectionWithReconstructionOperatorTraits>;
}; // class AdvectionWithReconstructionOperatorTraits


} // namespace internal


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class Traits =
              internal::AdvectionWithReconstructionOperatorTraits<AdvectionOperatorImp, ReconstructionOperatorImp>>
class AdvectionWithReconstructionOperator : public OperatorInterface<Traits>
{
public:
  using AdvectionOperatorType = typename Traits::AdvectionOperatorType;
  using ReconstructionOperatorType = typename Traits::ReconstructionOperatorType;

  AdvectionWithReconstructionOperator(const AdvectionOperatorType& advection_operator,
                                      const ReconstructionOperatorType& reconstruction_operator)
    : advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    using ReconstructedFunctionType = ReconstructedLocalizableFunction<typename SourceType::SpaceType::GridLayerType,
                                                                       typename SourceType::DomainFieldType,
                                                                       SourceType::dimDomain,
                                                                       typename SourceType::RangeFieldType,
                                                                       SourceType::dimRange,
                                                                       SourceType::dimRangeCols>;

    // do reconstruction
    ReconstructedFunctionType reconstructed_function(source.space().grid_layer());
    reconstruction_operator_.apply(source, reconstructed_function, param);

    // apply advection operator
    std::fill(range.vector().begin(), range.vector().end(), 0.);
    advection_operator_.apply(reconstructed_function, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
}; // class AdvectionWithReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_HH
