// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/ESV2007.hh"
#include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy")
      return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
      return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return {7.23e-02, 1.82e-02, 4.54e-03, 1.14e-03};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
      return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {4.48e-01, 2.07e-01, 9.91e-02, 4.85e-02};
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {1.37, 1.28, 1.23, 1.21};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {7.70e-01, 5.22e-01, 3.62e-01, 2.53e-01};
    else if (type
             == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {2.35e+00, 3.23e+00, 4.50e+00, 6.32e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex,
                                                                                                 conforming>,
                                                                                         double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
