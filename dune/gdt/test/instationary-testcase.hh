// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TEST_INSTATIONARY_TESTCASE_HH
#define DUNE_GDT_TEST_INSTATIONARY_TESTCASE_HH

#include <dune/xt/common/exceptions.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider/eoc.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \tparam ProblemType has to provide a type SolutionType which
 *         defines the type of the solution of the problem.
 * TODO: choose suitable SolutionType for Problems (provide Interface?)
 */
template <class GridImp, class ProblemImp>
class InstationaryTestCase : public XT::Grid::EOCGridProvider<GridImp>
{
  typedef XT::Grid::EOCGridProvider<GridImp> EocBaseType;

public:
  typedef ProblemImp ProblemType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::SolutionType SolutionType;

public:
  template <class... Args>
  InstationaryTestCase(const double divide_t_end_by_this, Args&&... args)
    : EocBaseType(std::forward<Args>(args)...)
    , divide_t_end_by_this_(divide_t_end_by_this)
    , zero_()
  {
  }

  virtual ~InstationaryTestCase() = default;

  virtual const ProblemType& problem() const = 0;

  virtual void print_header(std::ostream& out = std::cout) const
  {
    out << "+===============================================================+\n"
        << "|+=============================================================+|\n"
        << "||  This is a GDT::Tests:InstationaryTestCase, please provide ||\n"
        << "||  a meaningful message by implementing `print_header()`      ||\n"
        << "|+=============================================================+|\n"
        << "+===============================================================+" << std::endl;
  }

  virtual bool provides_exact_solution() const
  {
    return false;
  }

  virtual std::bitset<GridImp::dimension> periodic_directions() const
  {
    return std::bitset<GridImp::dimension>();
  }

  virtual double t_end() const
  {
    return problem().t_end() / divide_t_end_by_this_;
  }

  virtual const std::shared_ptr<const SolutionType> exact_solution() const
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
  const double divide_t_end_by_this_;
  const std::shared_ptr<const SolutionType> zero_;
}; // class InstationaryTestCase


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_TESTCASE_HH
