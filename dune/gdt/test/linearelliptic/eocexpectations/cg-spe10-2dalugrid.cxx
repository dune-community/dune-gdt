#include <config.h>
#if HAVE_ALUGRID

#include "cg-spe10-2dalugrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {3.44e-02, 1.01e-02};
  else if (type == "H1_semi")
    return {1.47e-01, 7.78e-02};
  else if (type == "energy")
    return {1.88e-01, 1.00e-01};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {5.08e-02, 1.54e-02};
  else if (type == "H1_semi")
    return {2.01e-01, 1.04e-01};
  else if (type == "energy")
    return {2.30e-01, 1.25e-01};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
} // ... results(...)

} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_ALUGRID
