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
        return {1.76e-01, 3.94e-02};
      else
        return {1.77e-01, 4.49e-02, 1.24e-02, 2.97e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {6.83e-01, 3.78e-01};
      else
        return {6.82e-01, 4.19e-01, 2.05e-01, 1.02e-01};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {4.35e-01, 2.66e-01};
      else
        return {4.29e-01, 2.63e-01, 1.84e-01, 1.44e-01};
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
