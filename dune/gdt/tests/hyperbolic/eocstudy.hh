// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_HYPERBOLIC_EOCSTUDY_HH
#define DUNE_GDT_TESTS_HYPERBOLIC_EOCSTUDY_HH

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>

#include "../nonstationary-eocstudy.hh"
#include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Tests {


template< class TestCaseImp, class DiscretizerImp >
class HyperbolicEocStudy
  : public NonStationaryEocStudy< TestCaseImp, DiscretizerImp >
{
  typedef NonStationaryEocStudy< TestCaseImp, DiscretizerImp > BaseType;
public:
  using typename BaseType::TestCaseType;
  using typename BaseType::Discretizer;
  using typename BaseType::DiscretizationType;
  using typename BaseType::GridViewType;
  using typename BaseType::VectorType;

  // a perfect forwarding ctor did not do the job here, since it was not able to match the std::initializer_list: {"L2"}
  HyperbolicEocStudy(TestCaseType& test_case,
                     const std::vector< std::string > only_these_norms = {},
                     const std::string visualize_prefix = "")
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
    //     template class HyperbolicEocExpectations< TestCasesType, Discretizer::type, GridViewType::dimension >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder and Discretizer::type,
    // if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class HyperbolicEocExpectations< ... >'
    // to each test source using these results!
    return HyperbolicEocExpectations< TestCaseType, Discretizer::type, GridViewType::dimension >::rate(type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker, see the explanation above in expected_rate()!
    return HyperbolicEocExpectations< TestCaseType, Discretizer::type, GridViewType::dimension >::results(this->test_case_, type);
  }

  virtual std::vector< std::string > available_norms() const override final
  {
    return {"L1"};
  }

  virtual double compute_norm(const GridViewType& grid_view,
                              const VectorType& function,
                              const std::string type) const override final
  {
    if (type == "L1") {
      double norm = 0;
      for (size_t ii = 0; ii < function.size(); ++ii) {
        double spatial_integral = 0;
        const auto it_end = grid_view.template end< 0 >();
        for (auto it = grid_view.template begin< 0 >(); it != it_end; ++it) {
          const auto& entity = *it;
          double value = std::abs(function[ii].second[grid_view.indexSet().index(entity)]);
          spatial_integral += value*entity.geometry().volume();
        }
        const double dt = (ii == function.size() - 1) ? function[ii].first - function[ii-1].first : function[ii+1].first - function[ii].first;
        norm += dt*spatial_integral;
      }
      return norm;
    }
    else {
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Wrong type `" << type << "` requested (see `available_norms()`!");
    return 0.0;
    }
  } // ... compute_norm(...)

  virtual std::vector< std::string > available_estimators() const override final
  {
    return {};
  }

  virtual double estimate(const VectorType& /*vector*/, const std::string /*type*/) const override final
  {
    std::abort();
  }

}; // class HyperbolicEocStudy


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_HYPERBOLIC_EOCSTUDY_HH
