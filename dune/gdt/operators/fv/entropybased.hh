// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
#define DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH

#include <dune/xt/grid/walker/functors.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "reconstructed_function.hh"

namespace Dune {
namespace GDT {


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp,
          class Traits>
class EntropyBasedMomentAdvectionOperator;


namespace internal {


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp>
class EntropyBasedMomentAdvectionOperatorTraits
{
public:
  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;
  using RealizabilityLimiterType = RealizabilityLimiterImp;
  using RegularizationOperatorType = RegularizationOperatorImp;

  using derived_type = EntropyBasedMomentAdvectionOperator<AdvectionOperatorImp,
                                                           ReconstructionOperatorImp,
                                                           RealizabilityLimiterImp,
                                                           RegularizationOperatorImp,
                                                           EntropyBasedMomentAdvectionOperatorTraits>;
}; // class EntropyBasedMomentAdvectionOperatorTraits


} // namespace internal


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp,
          class Traits = internal::EntropyBasedMomentAdvectionOperatorTraits<AdvectionOperatorImp,
                                                                             ReconstructionOperatorImp,
                                                                             RealizabilityLimiterImp,
                                                                             RegularizationOperatorImp>>
class EntropyBasedMomentAdvectionOperator : public OperatorInterface<Traits>
{
public:
  using AdvectionOperatorType = typename Traits::AdvectionOperatorType;
  using RegularizationOperatorType = typename Traits::RegularizationOperatorType;
  using ReconstructionOperatorType = typename Traits::ReconstructionOperatorType;
  using RealizabilityLimiterType = typename Traits::RealizabilityLimiterType;

  EntropyBasedMomentAdvectionOperator(const AdvectionOperatorType& advection_operator,
                                      const ReconstructionOperatorType reconstruction_operator,
                                      const RegularizationOperatorType& regularization_operator,
                                      const RealizabilityLimiterType& realizability_limiter)
    : advection_operator_(advection_operator)
    , reconstruction_operator_(reconstruction_operator)
    , regularization_operator_(regularization_operator)
    , realizability_limiter_(realizability_limiter)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    using ReconstructedFunctionType = ReconstructedLocalizableFunction<typename SourceType::GridLayerType,
                                                                       typename SourceType::DomainFieldType,
                                                                       SourceType::dimDomain,
                                                                       typename SourceType::RangeFieldType,
                                                                       SourceType::dimRange,
                                                                       SourceType::dimRangeCols>;
    // solve optimization problems and regularize if necessary
    regularization_operator_.apply(source, source, param);

    // do reconstruction
    ReconstructedFunctionType reconstructed_function(source.space().grid_layer());
    reconstruction_operator_.apply(source, reconstructed_function, param);

    // perform realizability limiting
    realizability_limiter_.apply(reconstructed_function, reconstructed_function, param);

    advection_operator_.apply(reconstructed_function, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  const RegularizationOperatorType& regularization_operator_;
  const RealizabilityLimiterType& realizability_limiter_;
}; // class EntropyBasedMomentAdvectionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
