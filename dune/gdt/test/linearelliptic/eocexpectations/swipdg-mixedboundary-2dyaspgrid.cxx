#include <config.h>

#include "swipdg-mixedboundary-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {7.38e-02, 1.84e-02};
    else
      return {7.56e-02, 2.06e-02, 5.56e-03, 1.36e-03};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {3.66e-01, 1.72e-01};
    else
      return {3.76e-01, 1.95e-01, 9.89e-02, 4.67e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {8.85e-03, 2.07e-03};
    else
      return {8.82e-03, 2.02e-03, 5.24e-04, 1.42e-04};
  } else if (type == "H1_semi" || type == "energy") {
    if (test_case.num_refinements() == 1)
      return {6.85e-02, 3.56e-02};
    else
      return {6.61e-02, 3.27e-02, 1.81e-02, 1.03e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
