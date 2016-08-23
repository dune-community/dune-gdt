// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/spe10.hh"
#include "../swipdg-estimator-expectations.hh"


namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


// polorder 1, conforming

std::vector<double>
results_LinearEllipticSwipdgEstimatorExpectationsSpe10Model1TestCaseALUGrid22simplexconformingdouble1swipdg1(
    const LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>& /*test_case*/,
    const std::string type)
{
  if (type == "energy")
    return {8.38e-01, 4.02e-01};
  else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_nonconformity_ESV2007_id())
    return {2.74e+00, 1.84e+00};
  else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
    return {2.26e-11, 4.39e-12};
  else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_diffusive_flux_ESV2007_id())
    return {1.22e+00, 7.62e-01};
  else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
    return {3.00e+00, 1.99e+00};
  else if (type == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_id())
    return {3.59e+00, 4.95e+00};
  else if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
    return {1.99e+00, 1.61e+00};
  else if (type
           == "efficiency_" + LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007_alternative_summation_id())
    return {2.38e+00, 4.01e+00};
  else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace internal


// polorder 1, conforming

#if HAVE_EIGEN

template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1,
                                                XT::LA::ChooseBackend::eigen_sparse, anything>
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {

    if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
#ifndef NDEBUG
      return {1.08e-12, 2.60e-12};
#else
      return {1.07e-12, 2.60e-12};
#endif
    else
      return internal::
          results_LinearEllipticSwipdgEstimatorExpectationsSpe10Model1TestCaseALUGrid22simplexconformingdouble1swipdg1(
              test_case, type);
  }
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                     conforming>,
                                                                                             double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg, 1,
                                                         XT::LA::ChooseBackend::eigen_sparse>;


#endif // HAVE_EIGEN
#if HAVE_DUNE_ISTL


template <bool anything>
class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                    double, 1>,
                                                LinearElliptic::ChooseDiscretizer::swipdg, 1,
                                                XT::LA::ChooseBackend::istl_sparse, anything>
    : public internal::LinearEllipticSwipdgEstimatorExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
#ifndef NDEBUG
    if (type == LinearElliptic::SwipdgFluxreconstrutionEstimators::local_residual_ESV2007_id())
      return {2.26e-11, 4.45e-12};
#endif
    return internal::
        results_LinearEllipticSwipdgEstimatorExpectationsSpe10Model1TestCaseALUGrid22simplexconformingdouble1swipdg1(
            test_case, type);
  }
}; // LinearEllipticSwipdgEstimatorExpectations


template class LinearEllipticSwipdgEstimatorExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex,
                                                                                                     conforming>,
                                                                                             double, 1>,
                                                         LinearElliptic::ChooseDiscretizer::swipdg, 1,
                                                         XT::LA::ChooseBackend::istl_sparse>;


#endif // HAVE_DUNE_ISTL

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
