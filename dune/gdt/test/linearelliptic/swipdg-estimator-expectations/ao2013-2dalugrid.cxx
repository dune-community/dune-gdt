// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/AO2013.hh"
#include "../swipdg-estimator-expectations.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg,
                                                1,
                                                anything>
{
  typedef LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
      return {1.58e-01, 2.79e-01};
#else
      return {9.10e-02, 5.23e-02, 2.68e-02, 1.20e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.02e-01, 2.19e-01};
#else
      return {9.57e-02, 1.10e-01, 5.12e-02, 2.17e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id()) {
      return {0.0, 0.0, 0.0, 0.0};
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.79e-01, 1.35e-01};
#else
      return {1.12e-01, 6.54e-02, 3.54e-02, 1.90e-02};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {1.48e-01, 1.28e-01, 6.22e-02, 2.89e-02};
#endif
    } else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {1.30e+00, 9.24e-01};
#else
      return {1.62e+00, 2.45e+00, 2.32e+00, 2.40e+00};
#endif
    } else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {};
#else
      return {4.56e-01, 4.19e-01, 2.94e-01, 2.02e-01};
#endif
    } else if (type
               == "efficiency_"
                      + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id()) {
#if DXT_DISABLE_LARGE_TESTS
      return {3.36e+00, 2.14e+00};
#else
      return {5.01e+00, 8.01e+00, 1.10e+01, 1.68e+01};
#endif
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  }
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::
                                                             AO2013TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                            double,
                                                                            1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg,
                                                         1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
