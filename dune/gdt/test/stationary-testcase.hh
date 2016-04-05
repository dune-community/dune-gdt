// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_STATIONARY_TESTCASE_HH
#define DUNE_GDT_TEST_STATIONARY_TESTCASE_HH

#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/provider/eoc.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \tparam ProblemType has to provide a type FunctionType (derived from Stuff::LocalizableFunctionInterface) which
 *         defines the type of the solution of the problem.
 */
template< class GridImp, class ProblemImp >
class StationaryTestCase
  : public Stuff::Grid::Providers::EOC< GridImp >
{
  typedef Stuff::Grid::Providers::EOC< GridImp > EocBaseType;
public:
  typedef ProblemImp                         ProblemType;
  typedef typename ProblemType::FunctionType FunctionType;
private:
  static_assert(Stuff::is_localizable_function< FunctionType >::value,
                "ProblemImp::FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef Stuff::Functions::Constant< typename FunctionType::EntityType,
                                      typename FunctionType::DomainFieldType, FunctionType::dimDomain,
                                      typename FunctionType::RangeFieldType, FunctionType::dimRange >
      ConstantFunctionType;

public:
  template< class... Args >
  StationaryTestCase(Args&& ...args)
    : EocBaseType(std::forward< Args >(args)...)
    , zero_(0.0)
  {}

  virtual ~StationaryTestCase() {}

  virtual const ProblemType& problem() const = 0;

  virtual void print_header(std::ostream& out = std::cout) const
  {
    out << "+==============================================================+\n"
        << "|+============================================================+|\n"
        << "||  This is a GDT::Tests::StationaryTestCase, please provide  ||\n"
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
      DUNE_THROW(Stuff::Exceptions::you_have_to_implement_this,
                 "If provides_exact_solution() is true, exact_solution() has to be implemented!");
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
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
