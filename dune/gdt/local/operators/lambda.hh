// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH
#define DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH

#include <functional>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards, needed for the Traits
template <class S, class RS, class RV>
class LocalLambdaOperator;

template <class S, class I, class RS, class RV>
class LocalLambdaCouplingOperator;


namespace internal {


template <class S, class RS, class RV>
class LocalLambdaOperatorTraits
{
  static_assert(XT::Functions::is_localizable_function<S>::value, "");
  static_assert(is_space<RS>::value, "");
  static_assert(XT::LA::is_vector<RV>::value, "");

public:
  typedef LocalLambdaOperator<S, RS, RV> derived_type;
};


template <class S, class I, class RS, class RV>
class LocalLambdaCouplingOperatorTraits
{
  static_assert(XT::Functions::is_localizable_function<S>::value, "");
  static_assert(XT::Grid::is_intersection<I>::value, "");
  static_assert(is_space<RS>::value, "");
  static_assert(XT::LA::is_vector<RV>::value, "");

public:
  using derived_type = LocalLambdaCouplingOperator<S, I, RS, RV>;
};


} // namespace internal


template <class S, class RS, class RV>
class LocalLambdaOperator : public LocalOperatorInterface<internal::LocalLambdaOperatorTraits<S, RS, RV>>
{
  typedef LocalLambdaOperator<S, RS, RV> ThisType;

public:
  typedef S SourceType;
  typedef LocalDiscreteFunction<RS, RV> LocalRangeType;
  typedef std::function<void(const SourceType&, LocalRangeType&)> LocalLambdaType;

  LocalLambdaOperator(LocalLambdaType local_lambda)
    : local_lambda_(local_lambda)
  {
  }

  template <class UnknownSourceType, class UnknownLocalRangeType>
  void apply(const UnknownSourceType& /*source*/, UnknownLocalRangeType& /*local_range*/) const
  {
    static_assert(AlwaysFalse<UnknownSourceType>::value, "");
    DUNE_THROW(NotImplemented, "For this unknown source and unknown local_range!");
  }

  void apply(const SourceType& source, LocalRangeType& local_range) const
  {
    local_lambda_(source, local_range);
  }

  LocalLambdaOperator(const ThisType&) = default;
  LocalLambdaOperator(ThisType&&) = default;

private:
  const LocalLambdaType local_lambda_;
}; // class LocalLambdaOperator


template <class S, class I, class RS, class RV>
class LocalLambdaCouplingOperator
    : public LocalCouplingOperatorInterface<internal::LocalLambdaCouplingOperatorTraits<S, I, RS, RV>>
{
  using ThisType = LocalLambdaCouplingOperator<S, I, RS, RV>;
  using BaseType = LocalCouplingOperatorInterface<internal::LocalLambdaCouplingOperatorTraits<S, I, RS, RV>>;

public:
  using Traits = internal::LocalLambdaCouplingOperatorTraits<S, I, RS, RV>;
  using SourceType = S;
  using IntersectionType = I;
  using LocalRangeType = LocalDiscreteFunction<RS, RV>;
  using LocalLambdaType =
      std::function<void(const SourceType&, const IntersectionType, LocalRangeType&, LocalRangeType&)>;

  LocalLambdaCouplingOperator(LocalLambdaType local_lambda)
    : local_lambda_(local_lambda)
  {
  }

  template <class UnknownSourceType, class UnknownIntersectionType, class UnknownRangeType>
  void apply(const UnknownSourceType& /*source*/,
             const UnknownIntersectionType& /*intersection*/,
             UnknownRangeType& /*local_range_entity*/,
             UnknownRangeType& /*local_range_neighbor*/) const
  {
    static_assert(AlwaysFalse<UnknownSourceType>::value, "");
    DUNE_THROW(NotImplemented, "For these unknown argument types!");
  }

  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalRangeType& local_range_entity,
             LocalRangeType& local_range_neighbor) const
  {
    local_lambda_(source, intersection, local_range_entity, local_range_neighbor);
  }

  LocalLambdaCouplingOperator(const ThisType&) = default;
  LocalLambdaCouplingOperator(ThisType&&) = default;

private:
  const LocalLambdaType local_lambda_;
}; // class LocalLambdaCouplingOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH
