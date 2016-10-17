#include <config.h>
#include "cg-mixedboundary-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {9.48e-02, 2.13e-02};
    else
      return {9.77e-02, 2.49e-02, 6.16e-03, 1.38e-03};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {3.64e-01, 1.67e-01};
    else
      return {3.71e-01, 1.87e-01, 9.31e-02, 4.31e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
