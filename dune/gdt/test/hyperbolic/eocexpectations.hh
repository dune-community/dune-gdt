// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2015 - 2016)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH

#include <vector>

#include "eocexpectations_base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class TestCaseType,
          Hyperbolic::ChooseDiscretizer disc,
          size_t dimDomain,
          NumericalFluxes num_flux,
          TimeStepperMethods time_stepper,
          TimeStepperMethods rhs_time_stepper>
class HyperbolicEocExpectations : public internal::HyperbolicEocExpectationsBase<dimDomain>
{
public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "  TestCaseType: " << XT::Common::Typename<TestCaseType>::value() << "\n"
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


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_EOCEXPECTATIONS_HH
