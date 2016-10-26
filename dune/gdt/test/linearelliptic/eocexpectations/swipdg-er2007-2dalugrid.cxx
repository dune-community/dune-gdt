#include <config.h>

#if HAVE_ALUGRID

#include "swipdg-er2007-2dalugrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


// polorder 1, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {1.11e-01, 6.09e-02, 1.65e-02};
  else if (type == "H1_semi" || type == "energy")
    return {3.66e-01, 2.98e-01, 1.46e-01};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, conforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {6.35e-02, 6.42e-03, 8.23e-04};
  else if (type == "H1_semi" || type == "energy")
    return {2.32e-01, 5.40e-02, 1.41e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 1, nonconforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 1>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {2.00e-01, 1.11e-01, 3.29e-02};
  else if (type == "H1_semi" || type == "energy")
    return {4.41e-01, 4.05e-01, 2.04e-01};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

// polorder 2, nonconforming

std::vector<double>
LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                              LinearElliptic::ChooseDiscretizer::swipdg, 2>::
    results(const LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 2>::TestCaseType&,
            const std::string type)
{
  if (type == "L2")
    return {1.51e-01, 1.27e-02, 1.50e-03};
  else if (type == "H1_semi" || type == "energy")
    return {3.69e-01, 8.28e-02, 2.12e-02};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
