// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Tobias Leibner  (2015)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include "discretizers/base.hh"

#include "eocexpectations_base.hh"

#include <dune/gdt/test/grids.hh>
#include "problems.hh"


namespace Dune {
namespace GDT {
namespace Test {

template <class TestCaseType, LinearElliptic::ChooseDiscretizer disc, int polOrder>
class LinearEllipticEocExpectations : public internal::LinearEllipticEocExpectationsBase<polOrder>
{
public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "  TestCaseType: " << XT::Common::Typename<TestCaseType>::value() << "\n"
                       << "  ChooseDiscretizer: " << int(disc) << "\n"
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


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_EOCEXPECTATIONS_HH
