#ifndef DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH
#define DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH

#include "discretizers/base.hh"

#include <dune/gdt/test/grids.hh>

namespace Dune {
namespace GDT {
namespace Test {


using Yasp1 = Yasp1Grid;
using Yasp2 = Yasp2Grid;

namespace internal {


template <int dimDomain>
class HyperbolicEocExpectationsBase
{
public:
  static double rate(const std::string type)
  {
    if (type == "L1") {
      if (dimDomain == 1)
        return 0.5;
      else
        return 0.25;
    } else {
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
      return 0;
    }
  } // ... rate(...)
}; // class HyperbolicEocExpectationsBase


} // namespace internal

template <class TestCaseType, Hyperbolic::ChooseDiscretizer disc, size_t dimDomain, NumericalFluxes num_flux,
          TimeStepperMethods time_stepper>
class HyperbolicEocExpectations;


} // namespace Test
} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH
