#include <config.h>
#if HAVE_ALUGRID

#include "cg-esv2007-2dalugrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {3.82e-02, 9.64e-03, 2.42e-03, 6.04e-04};
  else if (type == "H1_semi" || type == "energy")
    return {1.84e-01, 9.24e-02, 4.63e-02, 2.31e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::cg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::cg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {4.22e-02, 1.08e-02, 2.70e-03, 6.76e-04};
  else if (type == "H1_semi" || type == "energy")
    return {1.94e-01, 9.79e-02, 4.91e-02, 2.45e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
