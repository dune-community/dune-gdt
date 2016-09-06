#include <config.h>
#include "cg-spe10-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {1.86e-02, 1.51e-02};
  else if (type == "H1_semi")
    return {3.31e-01, 4.32e-01};
  else if (type == "energy")
    return {9.58e-01, 1.37e+00};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
