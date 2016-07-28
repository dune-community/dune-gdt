// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/dgf.hh>

#include "../problems/AO2013.hh"
#include "../swipdg-estimator-expectations.hh"

namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::AO2013TestCase<Alu2NonConformSimplex,
                                                                               double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<Alu2NonConformSimplex, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy")
      return {1.32e-01, 7.47e-02, 3.81e-02, 1.70e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
      return {1.30e-01, 1.42e-01, 6.30e-02, 2.65e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return {4.81e-16, 1.04e-14, 2.56e-15, 1.88e-13};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
      return {1.56e-01, 8.93e-02, 4.87e-02, 2.61e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {1.48e-01, 1.28e-01, 6.22e-02, 2.89e-02};
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {1.54e+00, 2.24e+00, 2.09e+00, 2.18e+00};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {4.56e-01, 4.19e-01, 2.94e-01, 2.02e-01};
    else if (type
             == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {4.05e+00, 6.44e+00, 8.77e+00, 1.35e+01};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::AO2013TestCase<Alu2NonConformSimplex,
                                                                                        double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
