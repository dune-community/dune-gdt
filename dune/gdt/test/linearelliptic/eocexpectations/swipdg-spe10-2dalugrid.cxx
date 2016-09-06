#include <config.h>

#if HAVE_ALUGRID

#include "swipdg-spe10-2dalugrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {9.48e-03, 2.64e-03};
  else if (type == "H1_semi")
    return {1.09e-01, 5.36e-02};
  else if (type == "energy")
    return {1.37e-01, 6.60e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {9.43e-03, 2.55e-03};
  else if (type == "H1_semi")
    return {1.43e-01, 7.04e-02};
  else if (type == "energy")
    return {1.74e-01, 8.53e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {2.74e-03, 8.75e-04};
  else if (type == "H1_semi")
    return {4.91e-02, 2.26e-02};
  else if (type == "energy")
    return {5.60e-02, 2.58e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {2.13e-03, 6.56e-04};
  else if (type == "H1_semi")
    return {6.41e-02, 3.24e-02};
  else if (type == "energy")
    return {7.56e-02, 3.81e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
#endif // HAVE_ALUGRID
