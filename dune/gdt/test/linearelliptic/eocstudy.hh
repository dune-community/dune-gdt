// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH

#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/gdt/operators/elliptic.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/operators/laplace.hh>

#include "../stationary-eocstudy.hh"
#include "eocexpectations.hh"

#include "eocexpectations/all.hh"


namespace Dune {
namespace GDT {
namespace Test {

template <class TestCaseImp, class DiscretizerImp>
class LinearEllipticEocStudy : public StationaryEocStudy<TestCaseImp, DiscretizerImp>
{
  typedef StationaryEocStudy<TestCaseImp, DiscretizerImp> BaseType;

public:
  using typename BaseType::TestCaseType;
  using typename BaseType::Discretizer;
  using typename BaseType::DiscretizationType;
  using typename BaseType::GridLayerType;
  using typename BaseType::FunctionType;
  using typename BaseType::VectorType;

private:
  static const int polOrder = Discretizer::polOrder;

public:
  // a perfect forwarding ctor did not do the job here, since it was not able to match the std::initializer_list: {"L2"}
  LinearEllipticEocStudy(TestCaseType& test_case,
                         const std::vector<std::string> only_these_norms = {},
                         const std::string visualize_prefix = "",
                         const size_t over_integrate = 2)
    : BaseType(test_case, only_these_norms, visualize_prefix)
    , over_integrate_(over_integrate)
  {
  }

  virtual ~LinearEllipticEocStudy() = default;

  virtual std::string identifier() const override
  {
    return Discretizer::static_id();
  }

  virtual size_t expected_rate(const std::string type) const override
  {
    // If you get an undefined reference here from the linker you are missing the appropriate
    // specialization of LinearEllipticEocExpectations!
    // For a new TestCaseType you have to add a specialization in a separate object file
    // (see linearelliptic-block-swipdg-expectations_os2014_2daluconform.cxx for example) and adjust the
    // CMakeLists.txt accordingly. For a new polOrder or Discretizer::type add
    //     template class LinearEllipticEocExpectations< TestCasesType, Discretizer::type, polOrder >;
    // in the appropriate (existing) object file and implement a specialization for this polOrder and Discretizer::type,
    // if needed!
    //
    // Oh: and do not forget to add
    //   'extern template class LinearEllipticEocExpectations< ... >'
    // to each test source using these results!
    return LinearEllipticEocExpectations<TestCaseType, Discretizer::type, polOrder>::rate(type);
  } // ... expected_rate(...)

  virtual std::vector<double> expected_results(const std::string type) const override
  {
    // If you get an undefined reference here from the linker, see the explanation above in expected_rate()!
    return LinearEllipticEocExpectations<TestCaseType, Discretizer::type, polOrder>::results(this->test_case_, type);
  }

  virtual std::vector<std::string> available_norms() const override
  {
    return {"L2", "H1_semi", "energy"};
  }

  virtual double
  compute_norm(const GridLayerType& grid_layer, const FunctionType& function, const std::string type) override final
  {
    if (type == "L2")
      return make_l2_operator(grid_layer, over_integrate_)->induced_norm(function);
    else if (type == "H1_semi")
      return make_laplace_operator(grid_layer, over_integrate_)->induced_norm(function);
    else if (type == "energy")
      return make_elliptic_operator(grid_layer,
                                    this->test_case_.problem().diffusion_factor(),
                                    this->test_case_.problem().diffusion_tensor(),
                                    over_integrate_)
          ->induced_norm(function);
    else
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                 "Wrong type `" << type << "` requested (see `available_norms()`!");
  } // ... compute_norm(...)

  virtual std::vector<std::string> available_estimators() const override
  {
    return {};
  }

  virtual double estimate(const VectorType& /*vector*/, const std::string /*type*/) override
  {
    DUNE_THROW(NotImplemented, "Do not call me if available_estimators().size() == 0!");
  }

private:
  const size_t over_integrate_;
}; // class LinearEllipticEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH
