#include <config.h>
#include "eocexpectations-fv-boltzmanncheckerboard-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Yasp2, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::Boltzmann2DCheckerboardTestCase<Yasp2, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 3.2))
      return {2.27e+00, 1.09e+00};
    else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 3.2 / 5.0))
      return {1.12e-01, 6.11e-02};
    else
      EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
