// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

# include <dune/grid/alugrid.hh>

# include "problems/AO2013.hh"
# include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {7.19e-02, 3.26e-02};
      else
        return {7.92e-02, 4.15e-02, 1.19e-02, 2.72e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {3.15e-01, 2.57e-01};
      else
        return {3.51e-01, 3.02e-01, 1.63e-01, 7.51e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {2.58e-01, 1.41e-01};
      else
        return {2.72e-01, 1.66e-01, 8.90e-02, 4.27e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {2.01e-01, 7.41e-02};
      else
        return {2.22e-01, 9.90e-02, 2.96e-02, 6.62e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {6.71e-01, 4.24e-01};
      else
        return {7.00e-01, 4.89e-01, 2.69e-01, 1.25e-01};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {4.03e-01, 2.13e-01};
      else
        return {4.21e-01, 2.50e-01, 1.34e-01, 6.36e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations



template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::AO2013TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
