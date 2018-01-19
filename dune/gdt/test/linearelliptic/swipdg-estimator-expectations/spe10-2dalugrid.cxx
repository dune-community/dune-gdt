// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#include "config.h"

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/grid.hh>

#include "../problems/spe10.hh"
#include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1,
                                                anything>
{
  typedef LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
      return {3.88e+01, 4.98e+01};
#else
      return {8.38e-01, 4.02e-01};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {3.20e+01, 5.47e+01};
#else
      return {2.74e+00, 1.84e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id()) {
      return {0.0, 0.0};
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {9.93e+00, 9.07e+00};
#else
      return {1.22e+00, 7.62e-01};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {3.00e+00, 1.99e+00};
#endif
    } else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {8.65e-01, 1.11e+00};
#else
      return {3.59e+00, 4.95e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {1.99e+00, 1.61e+00};
#endif
    } else if (type
               == "efficiency_"
                      + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.67e-01, 1.60e-01};
#else
      return {2.38e+00, 4.01e+00};
#endif
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  }
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                             Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg,
                                                         1>;


// polorder 1, nonconforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1,
                                                anything>
{
  typedef LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
      return {3.81e+01, 2.72e+01};
#else
      return {1.06e+00, 5.18e-01};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {3.20e+01, 4.28e+01};
#else
      return {2.90e+00, 1.91e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id()) {
      return {0.0, 0.0};
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {9.93e+00, 8.07e+00};
#else
      return {1.65e+00, 1.02e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {};
#endif
    } else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {8.79e-01, 1.60e+00};
#else
      return {3.16e+00, 4.18e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {};
#endif
    } else if (type
               == "efficiency_"
                      + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.70e-01, 2.62e-01};
#else
      return {2.02e+00, 3.30e+00};
#endif
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  }
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                             Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg,
                                                         1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
