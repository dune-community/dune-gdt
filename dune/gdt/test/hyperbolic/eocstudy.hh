// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_EOCSTUDY_HH
#define DUNE_GDT_TEST_HYPERBOLIC_EOCSTUDY_HH

#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include "../instationary-eocstudy.hh"
#include "eocexpectations_base.hh"
#include "all_eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class TestCaseImp, class DiscretizerImp>
class HyperbolicEocStudy : public InstationaryEocStudy<TestCaseImp, DiscretizerImp>
{
  typedef InstationaryEocStudy<TestCaseImp, DiscretizerImp> BaseType;

public:
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::DiscretizationType;
  using typename BaseType::Discretizer;
  using typename BaseType::GridLayerType;
  using typename BaseType::TestCaseType;

  // a perfect forwarding ctor did not do the job here, since it was not able to match the std::initializer_list: {"L2"}
  HyperbolicEocStudy(TestCaseType& test_case,
                     const std::vector<std::string> only_these_norms = {},
                     const std::string visualize_prefix = TestCaseType::ProblemType::static_id())
    : BaseType(test_case, only_these_norms, visualize_prefix)
  {}

  virtual ~HyperbolicEocStudy() {}

  virtual std::string identifier() const override final
  {
    return Discretizer::static_id();
  }

  virtual size_t expected_rate(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of HyperbolicEocExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder or Discretizer::type add
    //     template class HyperbolicEocExpectations< TestCasesType, Discretizer::type, GridLayerType::dimension >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder and Discretizer::type,
    // if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class HyperbolicEocExpectations< ... >'
    // to each test source using these results!
    return HyperbolicEocExpectations<TestCaseType,
                                     Discretizer::type,
                                     GridLayerType::dimension,
                                     Discretizer::numerical_flux_type,
                                     Discretizer::time_stepper_type,
                                     Discretizer::rhs_time_stepper_type,
                                     Discretizer::pol_order>::rate(type);
  } // ... expected_rate(...)

  virtual std::vector<double> expected_results(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker, see the explanation above in expected_rate()!
    return HyperbolicEocExpectations<TestCaseType,
                                     Discretizer::type,
                                     GridLayerType::dimension,
                                     Discretizer::numerical_flux_type,
                                     Discretizer::time_stepper_type,
                                     Discretizer::rhs_time_stepper_type,
                                     Discretizer::pol_order>::results(this->test_case_, type);
  }

  virtual std::vector<std::string> available_norms() const override final
  {
    return {"L1"};
  }

  virtual double compute_norm(const DiscreteSolutionType& solution, const std::string type) override final
  {
    if (type == "L1") {
      double norm = 0;
      // walk over time steps
      const auto solution_it_end = solution.end();
      for (auto solution_it = solution.begin(); solution_it != solution_it_end; ++solution_it) {
        double spatial_integral = 0;
        // walk over all entities, solution is constant on each entity
        const auto& grid_layer = solution_it->second.space().grid_layer();
        for (const auto& entity : Dune::elements(grid_layer)) {
          double value = 0;
          for (const auto& index : solution_it->second.space().mapper().globalIndices(entity))
            value += std::abs(solution_it->second.vector().get_entry(index));
          spatial_integral += value * entity.geometry().volume();
        }
        auto solution_it_copy = solution_it;
        const double dt = (solution_it == --solution.end()) ? ((*solution_it).first - (*(--solution_it_copy)).first)
                                                            : ((*(++solution_it_copy)).first - (*solution_it).first);
        norm += dt * spatial_integral;
      }
      return norm;
    } else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                 "Wrong type `" << type << "` requested (see `available_norms()`!");
      return 0.0;
    }
  } // ... compute_norm(...)

  virtual std::vector<std::string> available_estimators() const override final
  {
    return {};
  }

  virtual double estimate(const DiscreteSolutionType& /*solution*/, const std::string /*type*/) override final
  {
    std::abort();
  }

}; // class HyperbolicEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_HYPERBOLIC_EOCSTUDY_HH
