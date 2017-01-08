#include <config.h>
#include "eocexpectations-fv-sourcebeam-1dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 4.0))
      return {3.25e-01, 1.63e-01};
    else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 4.0 / 5.0))
      return {9.16e-02, 4.08e-02};
    else
      EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov_with_reconstruction,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::SourceBeamTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov_with_reconstruction,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 4.0))
        return {2.63e-01, 1.39e-01};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 4.0 / 5.0))
        return {7.07e-02, 3.13e-02};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
