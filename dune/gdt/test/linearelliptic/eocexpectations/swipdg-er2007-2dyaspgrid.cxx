#include <config.h>

#include "swipdg-er2007-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.19e+00, 2.30e-01};
#else
    return {2.30e-01, 1.60e-01, 4.88e-02};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.02e+00, 4.58e-01};
#else
    return {4.58e-01, 4.41e-01, 2.26e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.49e+00, 2.05e-01};
#else
    return {2.05e-01, 1.66e-02, 1.91e-03};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {7.98e-01, 4.39e-01};
#else
    return {4.39e-01, 9.32e-02, 2.37e-02};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
