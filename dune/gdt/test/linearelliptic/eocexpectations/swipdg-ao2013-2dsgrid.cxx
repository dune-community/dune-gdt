// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/sgrid.hh>

#include "../problems/AO2013.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {1.08e-01, 4.64e-02};
      else
        return {1.05e-01, 5.00e-02, 1.16e-02, 3.22e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {6.95e-01, 4.28e-01};
      else
        return {6.90e-01, 4.81e-01, 2.28e-01, 1.22e-01};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {5.22e-01, 3.57e-01};
      else
        return {5.09e-01, 3.44e-01, 2.64e-01, 2.20e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 2 >
{
  typedef LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {8.95e-02, 8.74e-03};
      else
        return {8.97e-02, 8.75e-03, 1.89e-03, 6.60e-04};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {5.28e-01, 1.51e-01};
      else
        return {5.24e-01, 1.47e-01, 7.39e-02, 5.43e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {2.91e-01, 2.20e-01};
      else
        return {2.70e-01, 2.05e-01, 1.77e-01, 1.46e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;


} // namespace Test
} // namespace GDT
} // namespace Dune
