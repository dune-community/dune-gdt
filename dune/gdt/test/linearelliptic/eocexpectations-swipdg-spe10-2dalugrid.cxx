// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

# include <dune/grid/alugrid.hh>

# include "problems/spe10.hh"
# include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {9.48e-03, 2.64e-03};
    else if (type == "H1_semi")
      return {1.09e-01, 5.36e-02};
    else if (type == "energy")
      return {1.37e-01, 6.60e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 2 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {};
    else if (type == "H1_semi")
      return {};
    else if (type == "energy")
      return {};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {3.19e-03, 9.01e-04};
    else if (type == "H1_semi")
      return {8.99e-02, 4.41e-02};
    else if (type == "energy")
      return {1.08e-01, 5.29e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::swipdg,
                                     2,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 2 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {};
    else if (type == "H1_semi")
      return {};
    else if (type == "energy")
      return {};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;
template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              2 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
