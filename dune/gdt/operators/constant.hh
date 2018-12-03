// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_CONSTANT_HH
#define DUNE_GDT_OPERATORS_CONSTANT_HH

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * See also OperatorInterface for a description of the template arguemnts.
 *
 * \sa OperatorInterface
 */
template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class ConstantOperator : public OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
{
  using ThisType = ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using BaseType = OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

public:
  using typename BaseType::F;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  ConstantOperator(const SourceSpaceType& src_space, const RangeSpaceType& rng_space, const VectorType& val)
    : source_space_(src_space)
    , range_space_(rng_space)
    , value_(val)
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
    DUNE_THROW_IF(!range_space_.contains(value_), Exceptions::operator_error, "");
    range = value_;
  }

  VectorType apply(const VectorType& source, const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!source_space_.contains(source), Exceptions::operator_error, "");
    DUNE_THROW_IF(!range_space_.contains(value_), Exceptions::operator_error, "");
    return value_;
  }

  std::vector<std::string> jacobian_options() const override final
  {
    return {"zero"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->jacobian_options().at(0), Exceptions::operator_error, "type = " << type);
    return {{"type", type}};
  }

  using BaseType::jacobian;

  void jacobian(const VectorType& /*source*/,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(
        jacobian_op.source_space().mapper().size() != source_space_.mapper().size(), Exceptions::operator_error, "");
    DUNE_THROW_IF(
        jacobian_op.range_space().mapper().size() != range_space_.mapper().size(), Exceptions::operator_error, "");
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, "opts = \n");
    DUNE_THROW_IF(
        opts.get<std::string>("type") != this->jacobian_options().at(0), Exceptions::operator_error, "opts = " << opts);
    // do nothing
  } // ... jacobian(...)

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const VectorType& value_;
}; // class ConstantOperator


template <class Matrix, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC, class V>
std::enable_if_t<XT::LA::is_matrix<Matrix>::value, ConstantOperator<Matrix, SGV, s_r, s_rC, r_r, r_rC, RGV>>
make_constant_operator(const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                       const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                       const XT::LA::VectorInterface<V>& value)
{
  return ConstantOperator<Matrix, SGV, s_r, s_rC, r_r, r_rC, RGV>(source_space, range_space, value.as_imp());
}


template <class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC, class V>
ConstantOperator<typename XT::LA::Container<F>::MatrixType, SGV, s_r, s_rC, r_r, r_rC, RGV>
make_constant_operator(const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                       const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                       const XT::LA::VectorInterface<V>& value)
{
  return make_constant_operator<typename XT::LA::Container<F>::MatrixType>(source_space, range_space, value.as_imp());
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_CONSTANT_HH
