// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include "../problems/ESV2007.hh"
#include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


// polorder 1

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<
    LinearElliptic::ESV2007TestCase<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1,
    anything> : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
      return {5.69e-1, 2.77e-1};
#else
      return {2.77e-01, 1.39e-01, 6.98e-02, 3.50e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.73e-1, 1.03e-1};
#else
      return {1.03e-01, 5.69e-02, 2.99e-02, 1.53e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {3.48e-1, 8.85e-2};
#else
      return {8.85e-02, 2.22e-02, 5.56e-03, 1.39e-03};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {8.97e-1, 4.42e-1};
#else
      return {4.42e-01, 2.23e-01, 1.12e-01, 5.64e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {2.16e-0, 1.91e-1};
#else
      return {0.0e-1, 0.0e-1};
#endif
    } else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {2.16e-0, 1.91e-0};
#else
      return {1.91e+00, 1.79e+00, 1.73e+00, 1.70e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {0.0e-1, 0.0e-1};
#else
      return {0.0e-1, 0.0e-1};
#endif
    } else if (type
               == "efficiency_"
                      + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {2.09e-0, 2.87e-0};
#else
      return {2.87e+00, 3.96e+00, 5.51e+00, 7.72e+00};
#endif
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {0.0e-1, 0.0e-1};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<
    LinearElliptic::ESV2007TestCase<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, double, 1>,
    LinearElliptic::ChooseDiscretizer::swipdg,
    1>;


} // namespace Test
} // namespace GDT
} // namespace Dune
