#include <config.h>
#include "cg-er2007-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {2.16e-01, 2.13e-01, 5.56e-02};
  else if (type == "H1_semi" || type == "energy")
    return {4.35e-01, 4.35e-01, 2.24e-01};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
