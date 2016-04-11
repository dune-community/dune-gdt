// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

# include <dune/grid/alugrid.hh>

# include "../problems/spe10.hh"
# include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template< bool anything >
class LinearEllipticSwipdgEstimatorExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                                 LinearElliptic::ChooseDiscretizer::swipdg,
                                                 1,
                                                 anything >
  : public internal::LinearEllipticSwipdgEstimatorExpectationsBase< 1 >
{
  typedef LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "energy")
      return {8.38e-01, 4.02e-01};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
      return {2.74e+00, 1.84e+00};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return {8.24e-12, 1.47e-11};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
      return {1.22e+00, 7.62e-01};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {3.00e+00, 1.99e+00};
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
      return {3.59e+00, 4.95e+00};
    else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {1.99e+00, 1.61e+00};
    else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
      return {2.38e+00, 4.01e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations< LinearElliptic::Spe10Model1TestCase< ALUGrid< 2, 2, simplex, conforming >, double, 1 >,
                                              LinearElliptic::ChooseDiscretizer::swipdg,
                                              1 >;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
