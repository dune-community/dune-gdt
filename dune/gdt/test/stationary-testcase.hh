// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)

#ifndef DUNE_GDT_TEST_STATIONARY_TESTCASE_HH
#define DUNE_GDT_TEST_STATIONARY_TESTCASE_HH

#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/gridprovider/eoc.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \tparam ProblemType has to provide a type FunctionType (derived from XT::Functions::LocalizableFunctionInterface)
 * which
 *         defines the type of the solution of the problem.
 */
template <class GridImp, class ProblemImp, typename DdGridImp = int>
class StationaryTestCase : public XT::Grid::EOCGridProvider<GridImp, DdGridImp>
{
  typedef XT::Grid::EOCGridProvider<GridImp, DdGridImp> EocBaseType;

public:
  typedef ProblemImp ProblemType;
  typedef typename ProblemType::FunctionType FunctionType;

private:
  static_assert(XT::Functions::is_localizable_function<FunctionType>::value,
                "ProblemImp::FunctionType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  typedef XT::Functions::ConstantFunction<typename FunctionType::EntityType,
                                          typename FunctionType::DomainFieldType,
                                          FunctionType::dimDomain,
                                          typename FunctionType::RangeFieldType,
                                          FunctionType::dimRange>
      ConstantFunctionType;

public:
  template <class... Args>
  StationaryTestCase(Args&&... args)
    : EocBaseType(std::forward<Args>(args)...)
    , zero_(0.0)
  {
  }

  virtual ~StationaryTestCase() = default;

  virtual const ProblemType& problem() const = 0;

  virtual void print_header(std::ostream& out = DXTC_LOG_INFO_0) const
  {
    out << "+==============================================================+\n"
        << "|+============================================================+|\n"
        << "||  This is a GDT::Test::StationaryTestCase, please provide   ||\n"
        << "||  a meaningful message by implementing `print_header()`     ||\n"
        << "|+============================================================+|\n"
        << "+==============================================================+" << std::endl;
  }

  virtual bool provides_exact_solution() const
  {
    return false;
  }

  virtual const FunctionType& exact_solution() const
  {
    if (provides_exact_solution())
      DUNE_THROW(XT::Common::Exceptions::you_have_to_implement_this,
                 "If provides_exact_solution() is true, exact_solution() has to be implemented!");
    else
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Do not call exact_solution() if provides_exact_solution() is false!");
    return zero_;
  }

private:
  const ConstantFunctionType zero_;
}; // class StationaryTestCase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_TESTCASE_HH
