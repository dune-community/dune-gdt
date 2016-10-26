#include <config.h>

#include "swipdg-esv2007-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {1.13e-02, 2.90e-03, 7.41e-04, 1.88e-04};
  else if (type == "H1_semi" || type == "energy")
    return {1.25e-01, 6.25e-02, 3.14e-02, 1.58e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {9.97e-04, 2.59e-04, 6.72e-05, 1.72e-05};
  else if (type == "H1_semi" || type == "energy")
    return {1.93e-02, 1.02e-02, 5.30e-03, 2.71e-03};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
