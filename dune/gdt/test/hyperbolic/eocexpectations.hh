// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include "discretizers/base.hh"

namespace Dune {
namespace GDT {
namespace Tests {
namespace internal {


template< int dimDomain >
class HyperbolicEocExpectationsBase
{
public:
  static double rate(const std::string type)
  {
    if (type == "L1") {
      if (dimDomain == 1)
        return 0.5;
      else
        return 0.25;
    } else {
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
      return 0;
    }
  } // ... rate(...)
}; // class HyperbolicEocExpectationsBase


} // namespace internal

template< class TestCaseType, Hyperbolic::ChooseDiscretizer disc, size_t dimDomain, Hyperbolic::FluxTimeStepperKombinations flux_timestepper, bool anything = true >
class HyperbolicEocExpectations
  : public internal::HyperbolicEocExpectationsBase< dimDomain >
{
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "  TestCaseType: " << Stuff::Common::Typename< TestCaseType >::value() << "\n"
                       << "  ChooseDiscretizer: ??\n"
                       << "  type: " << type << "\n"
                       << "Please put an appropriate specialization of HyperbolicEocExpectations for these\n"
                       << "combinations in a separate object file or add\n"
                       << "  'template class HyperbolicEocExpectations< ... >'\n"
                       << "for this ChooseDiscretizer in the appropriate object file!\n\n"
                       << "Oh: and do not forget to add\n"
                       << "  'extern template class HyperbolicEocExpectations< ... >'\n"
                       << "to each test source using these results!";
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH
