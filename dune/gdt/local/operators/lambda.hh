// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH
#define DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH

#include <functional>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, needed for the Traits
template <class S, class RS, class RV>
class LocalLambdaOperator;


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
    DUNE_THROW(NotImplemented, "For this unknown source and unknown local_range!");
  }

  void apply(const SourceType& source, LocalRangeType& local_range) const
  {
    local_lambda_(source, local_range);
  }

  LocalLambdaOperator(const ThisType&) = default;
  LocalLambdaOperator(ThisType&&) = default;

  const LocalLambdaType local_lambda_;
}; // class LocalLambdaOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_LAMBDA_HH
