// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_GENERIC_HH
#define DUNE_GDT_LOCAL_OPERATORS_GENERIC_HH

#include <functional>
#include <memory>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \note See LocalElementOperatorInterface for a description of the template arguments.
 *
 * \sa LocalElementOperatorInterface
 * \sa std::function
 */
template <class SV,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = SF,
          class RGV = SGV,
          class RV = SV>
class GenericLocalElementOperator : public LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>
{
  using ThisType = GenericLocalElementOperator;
  using BaseType = LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>;

public:
  using typename BaseType::LocalRangeType;
  using typename BaseType::LocalSourceType;
  using typename BaseType::SourceType;

  using GenericFunctionType = std::function<void(const SourceType& /*source*/,
                                                 const std::vector<std::unique_ptr<LocalSourceType>>& /*local_source*/,
                                                 LocalRangeType& /*local_range*/,
                                                 const XT::Common::Parameter& /*param*/)>;

  // When using this constructor, source has to be set by a call to with_source before calling apply
  GenericLocalElementOperator(GenericFunctionType func,
                              const size_t num_local_sources = 0,
                              const XT::Common::ParameterType& param_type = {})
    : BaseType(num_local_sources, param_type)
    , func_(func)
  {}

  GenericLocalElementOperator(const SourceType& source,
                              GenericFunctionType func,
                              const size_t num_local_sources = 0,
                              const XT::Common::ParameterType& param_type = {})
    : BaseType(source, num_local_sources, param_type)
    , func_(func)
  {}

  GenericLocalElementOperator(const ThisType& other)
    : BaseType(other)
    , func_(other.func_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const override final
  {
    func_(this->source(), this->local_sources(), local_range, this->parse_parameter(param));
  }

private:
  const GenericFunctionType func_;
}; // class GenericLocalElementOperator


/**
 * \note See LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 * \sa std::function
 */
template <class I,
          class SV,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = SF,
          class IRGV = SGV,
          class IRV = SV,
          class ORGV = IRGV,
          class ORV = IRV>
class GenericLocalIntersectionOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, IRGV, IRV, ORGV, ORV>
{
  using ThisType = GenericLocalIntersectionOperator;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, IRGV, IRV, ORGV, ORV>;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::LocalSourceType;
  using typename BaseType::SourceType;

  using GenericFunctionType = std::function<void(const SourceType& /*source*/,
                                                 const std::vector<std::unique_ptr<LocalSourceType>>& /*local_source*/,
                                                 const IntersectionType& /*intersection*/,
                                                 LocalInsideRangeType& /*local_range_inside*/,
                                                 LocalOutsideRangeType& /*local_range_outside*/,
                                                 const XT::Common::Parameter& /*param*/)>;

  // When using this constructor, source has to be set by a call to with_source before calling apply
  GenericLocalIntersectionOperator(GenericFunctionType func,
                                   const size_t num_local_sources = 1,
                                   const XT::Common::ParameterType& param_type = {})
    : BaseType(num_local_sources, param_type)
    , func_(func)
  {}

  GenericLocalIntersectionOperator(const SourceType& source,
                                   GenericFunctionType func,
                                   const size_t num_local_sources = 1,
                                   const XT::Common::ParameterType& param_type = {})
    : BaseType(source, num_local_sources, param_type)
    , func_(func)
  {}

  GenericLocalIntersectionOperator(const ThisType& other)
    : BaseType(other)
    , func_(other.func_)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& local_range_outside,
             const XT::Common::Parameter& param = {}) const override final
  {
    func_(this->source(),
          this->local_sources(),
          intersection,
          local_range_inside,
          local_range_outside,
          this->parse_parameter(param));
  }

private:
  const GenericFunctionType func_;
}; // class GenericLocalIntersectionOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_GENERIC_HH
