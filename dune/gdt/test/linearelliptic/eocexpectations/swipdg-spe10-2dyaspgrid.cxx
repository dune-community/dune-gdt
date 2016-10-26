#include <config.h>

#include "swipdg-spe10-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {2.73e-02, 2.35e-02};
  else if (type == "H1_semi")
    return {4.64e-01, 6.31e-01};
  else if (type == "energy")
    return {1.03e+00, 1.49e+00};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {1.82e-02, 7.50e-03};
  else if (type == "H1_semi")
    return {5.53e-01, 4.55e-01};
  else if (type == "energy")
    return {1.62e+00, 1.33e+00};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
