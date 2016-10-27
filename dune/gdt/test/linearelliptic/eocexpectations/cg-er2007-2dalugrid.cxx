#include <config.h>

#if HAVE_ALUGRID

#include "cg-er2007-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.74e+00, 2.16e-01};
#else
    return {2.14e-01, 1.48e-01, 3.82e-02};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.01e+00, 4.35e-01};
#else
    return {4.35e-01, 3.59e-01, 1.84e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.74e+00, 2.15e-01};
#else
    return {2.15e-01, 2.13e-01, 5.56e-02};
#endif
  } else if (type == "H1_semi" || type == "energy") {
#if DXT_DISABLE_LARGE_TESTS
    return {1.01e+00, 4.35e-01};
#else
    return {4.35e-01, 4.35e-01, 2.24e-01};
#endif
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_ALUGRID
