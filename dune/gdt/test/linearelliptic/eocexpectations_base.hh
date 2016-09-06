#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH

#include "discretizers/base.hh"
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/gdt/test/grids.hh>

namespace Dune {
namespace GDT {
namespace Test {

static const auto CG = LinearElliptic::ChooseDiscretizer::cg;

namespace internal {


template <int polOrder>
class LinearEllipticEocExpectationsBase
{
public:
  static size_t rate(const std::string type)
  {
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type == "energy")
      return polOrder;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class LinearEllipticEocExpectationsBase


} // namespace internal

template <class TestCaseType, LinearElliptic::ChooseDiscretizer disc, int polOrder>
class LinearEllipticEocExpectations;

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH
