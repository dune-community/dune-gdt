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


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
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
        return {4.39e-02, 1.11e-02};
      else
        return {4.52e-02, 1.27e-02, 3.28e-03, 8.50e-04};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {4.07e-01, 1.81e-01};
      else
        return {4.19e-01, 2.07e-01, 1.05e-01, 5.21e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {2.43e-01, 1.43e-01};
      else
        return {2.50e-01, 1.64e-01, 1.15e-01, 7.81e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< SGrid< 2, 2 >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune
