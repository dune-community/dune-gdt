// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH

#include <dune/stuff/test/gtest/gtest.h>

#include "discretizers/base.hh"

namespace Dune {
namespace GDT {
namespace Tests {
namespace internal {


template< int polOrder >
class LinearEllipticEocExpectationsBase
{
public:
  static size_t rate(const std::string type)
  {
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type == "energy")
      return polOrder;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class LinearEllipticEocExpectationsBase


} // namespace internal


template< class TestCaseType, LinearElliptic::ChooseDiscretizer disc, int polOrder, bool anything = true >
class LinearEllipticEocExpectations
  : public internal::LinearEllipticEocExpectationsBase< polOrder >
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "  TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
                       << "  ChooseDiscretizer: ??\n"
                       << "  polOrder: " << polOrder << "\n"
                       << "  type: " << type << "\n"
                       << "Please put an appropriate specialiaztion of LinearEllipticEocExpectations for these\n"
                       << "combinations in a separate object file or add\n"
                       << "  'template class LinearEllipticEocExpectations< ... >'\n"
                       << "for this ChooseDiscretizer and polOrder in the appropriate object file!\n\n"
                       << "Oh: and do not forget to add\n"
                       << "  'extern template class LinearEllipticEocExpectations< ... >'\n"
                       << "to each test source using these results!";
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH
