// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

//#if HAVE_ALUGRID // <- this is a tricky thing, since HAVE_ALUGRID is not defined here. This is the case since we cannot
                   //    add_dune_alugrid_flags(...) for this object file

#include <dune/grid/alugrid/common/declaration.hh>

#include "problems/ESV2007.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {4.23e-02, 9.65e-03, 2.42e-03, 6.05e-04};
    else if (type == "H1_semi" || type == "energy")
      return {4.32e-01, 2.05e-01, 1.03e-01, 5.14e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template< bool anything >
class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                     LinearElliptic::ChooseDiscretizer::cg,
                                     1,
                                     anything >
  : public internal::LinearEllipticEocExpectationsBase< 1 >
{
  typedef LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {4.23e-02, 1.08e-02, 2.70e-03, 6.76e-04};
    else if (type == "H1_semi" || type == "energy")
      return {4.32e-01, 2.18e-01, 1.09e-01, 5.45e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations



template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;
template class LinearEllipticEocExpectations< LinearElliptic::ESV2007TestCase< ALUGrid< 2, 2, simplex, nonconforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::cg,
                                              1 >;


} // namespace Tests
} // namespace GDT
} // namespace Dune

//#endif // HAVE_ALUGRID
