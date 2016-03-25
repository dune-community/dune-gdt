// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "problems/AO2013.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {4.92e-02, 1.09e-02};
      else
        return {5.02e-02, 1.18e-02, 3.37e-03, 1.06e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {4.71e-01, 2.09e-01};
      else
        return {4.83e-01, 2.32e-01, 1.28e-01, 7.95e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {3.53e-01, 2.77e-01};
      else
        return {3.36e-01, 2.53e-01, 2.06e-01, 1.74e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {};
      else
        return {};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<SGrid<2, 2>, double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune
