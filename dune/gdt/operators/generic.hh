// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_GENERIC_HH
#define DUNE_GDT_OPERATORS_GENERIC_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


/**
 * See also OperatorInterface for a description of the template arguemnts.
 *
 * \sa OperatorInterface
 */
template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class GenericOperator : public OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
{
  using ThisType = GenericOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using BaseType = OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

public:
  using typename BaseType::ConstSourceFunctionType;
  using typename BaseType::F;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;
  using LocalizableOperatorType =
      LocalizableOperatorBase<SGV, VectorType, s_r, s_rC, F, SGV, r_r, r_rC, F, RGV, VectorType>;
  using GenericElementOperatorType = GenericLocalElementOperator<VectorType, SGV, s_r, s_rC, F, r_r, r_rC>;
  using GenericIntersectionOperatorType =
      GenericLocalIntersectionOperator<XT::Grid::extract_intersection_t<SGV>, VectorType, SGV, s_r, s_rC, F, r_r, r_rC>;
  using GenericElementFunctionType = typename GenericElementOperatorType::GenericFunctionType;
  using GenericIntersectionFunctionType = typename GenericIntersectionOperatorType::GenericFunctionType;

  GenericOperator(const SourceSpaceType& src_space,
                  const RangeSpaceType& rng_space,
                  std::vector<GenericElementFunctionType> element_functions = {},
                  std::vector<GenericIntersectionFunctionType> intersection_functions = {})
    : source_space_(src_space)
    , range_space_(rng_space)
    , element_functions_(element_functions)
    , intersection_functions_(intersection_functions)
  {}

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  bool linear() const override final
  {
    return false;
  }

  using BaseType::apply;

  void
  apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!source_space_.contains(source), Exceptions::operator_error, "");
    DUNE_THROW_IF(!range_space_.contains(range), Exceptions::operator_error, "");
    ConstSourceFunctionType source_func = make_discrete_function(source_space_, source);
    RangeFunctionType range_func = make_discrete_function(range_space_, range);
    LocalizableOperatorType localizable_op(source_space().grid_view(), source_func, range_func);
    for (auto&& element_func : element_functions_)
      localizable_op.append(element_func);
    for (auto&& intersection_func : intersection_functions_)
      localizable_op.append(intersection_func);
    localizable_op.assemble(true);
  }

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::vector<GenericElementFunctionType> element_functions_;
  std::vector<GenericIntersectionFunctionType> intersection_functions_;
}; // class GenericOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_GENERIC_HH
