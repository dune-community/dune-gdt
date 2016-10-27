// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Tobias Leibner  (2015)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATOR_EXPECTATIONS_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATOR_EXPECTATIONS_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/container-interface.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include "discretizers/base.hh"
#include "estimators/swipdg-fluxreconstruction.hh"

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


template <int polOrder>
class LinearEllipticSwipdgEstimatorExpectationsBase
{
public:
  static size_t rate(const std::string type)
  {
    if (type == "energy")
      return polOrder;
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
      return polOrder;
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return polOrder + 1;
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
      return polOrder;
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return polOrder;
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return 0;
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return polOrder;
    else if (type
             == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return 0;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class LinearEllipticSwipdgEstimatorExpectationsBase


} // namespace internal


template <class TestCaseType, LinearElliptic::ChooseDiscretizer disc, int polOrder, bool anything = true>
class LinearEllipticSwipdgEstimatorExpectations
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<polOrder>
{
public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    EXPECT_TRUE(false) << "Please record the expected results for\n"
                       << "  TestCaseType: " << XT::Common::Typename<TestCaseType>::value() << "\n"
                       << "  ChooseDiscretizer: " << int(disc) << "\n"
                       << "  polOrder: " << polOrder << "\n"
                       << "  type: " << type << "\n"
                       << "Please put an appropriate specialiaztion of LinearEllipticSwipdgEstimatorExpectations for\n"
                       << "these combinations in a separate object file or add\n"
                       << "  'template class LinearEllipticSwipdgEstimatorExpectations< ... >'\n"
                       << "for this ChooseDiscretizer and polOrder in the appropriate object file!\n\n"
                       << "Oh: and do not forget to add\n"
                       << "  'extern template class LinearEllipticSwipdgEstimatorExpectations< ... >'\n"
                       << "to each test source using these results!";
    return {};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATOR_EXPECTATIONS_HH
