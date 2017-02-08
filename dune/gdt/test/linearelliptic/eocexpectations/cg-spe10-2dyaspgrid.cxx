#include <config.h>

#include "cg-spe10-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.41e+00, 1.00e+00};
#else
    return {1.86e-02, 1.51e-02};
#endif
  } else if (type == "H1_semi") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.10e+00, 1.00e+00};
#else
    return {3.31e-01, 4.32e-01};
#endif
  } else if (type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.28e+00, 1.00e+00};
#else
    return {9.58e-01, 1.37e+00};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
