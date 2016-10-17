#include <config.h>
#include "cg-ao2013-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L2") {
    if (test_case.num_refinements() == 1)
      return {1.76e-01, 3.94e-02};
    else
      return {1.77e-01, 4.49e-02, 1.24e-02, 2.97e-03};
  } else if (type == "H1_semi") {
    if (test_case.num_refinements() == 1)
      return {6.83e-01, 3.78e-01};
    else
      return {6.82e-01, 4.19e-01, 2.05e-01, 1.02e-01};
  } else if (type == "energy") {
    if (test_case.num_refinements() == 1)
      return {4.35e-01, 2.66e-01};
    else
      return {4.29e-01, 2.63e-01, 1.84e-01, 1.44e-01};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
