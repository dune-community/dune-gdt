// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/dgf.hh>

#include "../problems/ESV2007.hh"
#include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<Alu2NonConformSimplex,
                                                                                double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<Alu2NonConformSimplex, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy")
      return {3.33e-01, 1.63e-01, 8.07e-02, 4.01e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
      return {1.99e-01, 9.86e-02, 4.91e-02, 2.46e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return {7.23e-02, 1.82e-02, 4.54e-03, 1.14e-03};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
      return {3.47e-01, 1.71e-01, 8.45e-02, 4.20e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {4.48e-01, 2.07e-01, 9.91e-02, 4.85e-02};
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {1.37, 1.29, 1.25, 1.23};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {7.70e-01, 5.22e-01, 3.62e-01, 2.53e-01};
    else if (type
             == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {2.37e+00, 3.29e+00, 4.61e+00, 6.48e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<Alu2NonConformSimplex,
                                                                                         double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
