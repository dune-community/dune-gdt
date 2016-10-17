#include <config.h>
#if HAVE_ALUGRID

#include "cg-mixedboundary-2dalugrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {7.95e-02, 1.81e-02};
    else
      return {8.31e-02, 2.22e-02, 5.52e-03, 1.19e-03};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {3.01e-01, 1.42e-01};
    else
      return {3.11e-01, 1.64e-01, 8.23e-02, 3.75e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {1.28e-01, 3.14e-02};
    else
      return {1.35e-01, 3.99e-02, 1.02e-02, 2.16e-03};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {3.61e-01, 1.80e-01};
    else
      return {3.75e-01, 2.08e-01, 1.06e-01, 4.87e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_ALUGRID
