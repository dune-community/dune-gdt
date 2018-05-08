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

#include <dune/xt/grid/walker/functors.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "entropybased/realizability.hh"
#include "entropybased/regularization.hh"
#include "reconstruction/reconstructed_function.hh"

namespace Dune {
namespace GDT {


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp,
          class Traits>
class EntropyBasedMomentFvOperator;


namespace internal {


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp>
class EntropyBasedMomentFvOperatorTraits
{
public:
  using AdvectionOperatorType = AdvectionOperatorImp;
  using ReconstructionOperatorType = ReconstructionOperatorImp;
  using RealizabilityLimiterType = RealizabilityLimiterImp;
  using RegularizationOperatorType = RegularizationOperatorImp;
  using FieldType = typename AdvectionOperatorType::DomainFieldType;
  using JacobianType = NoJacobian;

  using derived_type = EntropyBasedMomentFvOperator<AdvectionOperatorImp,
                                                    ReconstructionOperatorImp,
                                                    RealizabilityLimiterImp,
                                                    RegularizationOperatorImp,
                                                    EntropyBasedMomentFvOperatorTraits>;
}; // class EntropyBasedMomentFvOperatorTraits


} // namespace internal


template <class AdvectionOperatorImp,
          class ReconstructionOperatorImp,
          class RealizabilityLimiterImp,
          class RegularizationOperatorImp,
          class Traits = internal::EntropyBasedMomentFvOperatorTraits<AdvectionOperatorImp,
                                                                      ReconstructionOperatorImp,
                                                                      RealizabilityLimiterImp,
                                                                      RegularizationOperatorImp>>
class EntropyBasedMomentFvOperator : public OperatorInterface<Traits>
{
public:
  using AdvectionOperatorType = typename Traits::AdvectionOperatorType;
  using RegularizationOperatorType = typename Traits::RegularizationOperatorType;
  using ReconstructionOperatorType = typename Traits::ReconstructionOperatorType;
  using RealizabilityLimiterType = typename Traits::RealizabilityLimiterType;

  EntropyBasedMomentFvOperator(const AdvectionOperatorType& advection_operator,
                               const ReconstructionOperatorType& reconstruction_operator,
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
    static_assert(is_discrete_function<SourceType>::value || is_const_discrete_function<SourceType>::value,
                  "SourceType has to be derived from (Const)DiscreteFunction!");
    static_assert(is_discrete_function<RangeType>::value, "RangeType has to be derived from DiscreteFunction!");
    using ReconstructedFunctionType = ReconstructedLocalizableFunction<typename SourceType::SpaceType::GridLayerType,
                                                                       typename SourceType::DomainFieldType,
                                                                       SourceType::dimDomain,
                                                                       typename SourceType::RangeFieldType,
                                                                       SourceType::dimRange,
                                                                       SourceType::dimRangeCols>;
    // solve optimization problems and regularize if necessary
    regularization_operator_.apply(source, range, param);

    // do reconstruction
    ReconstructedFunctionType reconstructed_function(range.space().grid_layer());
    reconstruction_operator_.apply(range, reconstructed_function, param);

    // perform realizability limiting
    realizability_limiter_.apply(range, reconstructed_function, param);

    std::fill(range.vector().begin(), range.vector().end(), 0.);
    advection_operator_.apply(reconstructed_function, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const ReconstructionOperatorType& reconstruction_operator_;
  const RegularizationOperatorType& regularization_operator_;
  const RealizabilityLimiterType& realizability_limiter_;
}; // class EntropyBasedMomentFvOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
