#include <config.h>

#include "cg-esv2007-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                  LinearElliptic::ChooseDiscretizer::cg,
                                                  1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg,
                                                1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {3.53e-02, 8.28e-03};
#else
    return {8.28e-03, 2.04e-03, 5.09e-04, 1.27e-04};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {2.33e-01, 1.14e-01};
#else
    return {1.14e-01, 5.68e-02, 2.83e-02, 1.42e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
