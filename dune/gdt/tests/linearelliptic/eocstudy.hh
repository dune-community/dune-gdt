// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>

#include "../stationary-eocstudy.hh"
#include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Tests {


template< class TestCaseImp, class DiscretizerImp >
class LinearEllipticEocStudy
  : public StationaryEocStudy< TestCaseImp, DiscretizerImp >
{
  typedef StationaryEocStudy< TestCaseImp, DiscretizerImp > BaseType;
public:
  using typename BaseType::TestCaseType;
  using typename BaseType::Discretizer;
  using typename BaseType::DiscretizationType;
  using typename BaseType::GridViewType;
  using typename BaseType::FunctionType;
  using typename BaseType::VectorType;
private:
  static const int polOrder = Discretizer::polOrder;
public:

  // a perfect forwarding ctor did not do the job here, since it was not able to match the std::initializer_list: {"L2"}
  LinearEllipticEocStudy(TestCaseType& test_case,
                         const std::vector< std::string > only_these_norms = {},
                         const std::string visualize_prefix = "",
                         const size_t over_integrate = 0)
    : BaseType(test_case, only_these_norms, visualize_prefix)
    , over_integrate_(over_integrate)
  {}

  virtual ~LinearEllipticEocStudy() {}

  virtual std::string identifier() const override final
  {
    return Discretizer::static_id();
  }

  virtual size_t expected_rate(const std::string type) const override final
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
    return LinearEllipticEocExpectations< TestCaseType, Discretizer::type, polOrder >::rate(type);
  } // ... expected_rate(...)

  virtual std::vector< double > expected_results(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker, see the explanation above in expected_rate()!
    return LinearEllipticEocExpectations< TestCaseType, Discretizer::type, polOrder >::results(this->test_case_, type);
  }

  virtual std::vector< std::string > available_norms() const override final
  {
    return {"L2", "H1_semi", "energy"};
  }

  virtual double compute_norm(const GridViewType& grid_view,
                              const FunctionType& function,
                              const std::string type) const override final
  {
    if (type == "L2")
      return Products::L2< GridViewType, double >(grid_view, over_integrate_).induced_norm(function);
    else if (type == "H1_semi")
      return Products::H1Semi< GridViewType, double >(grid_view, over_integrate_).induced_norm(function);
    else if (type == "energy")
      return Products::make_elliptic< double >(grid_view,
                                               this->test_case_.problem().diffusion_factor(),
                                               this->test_case_.problem().diffusion_tensor(),
                                               over_integrate_).induced_norm(function);
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Wrong type `" << type << "` requested (see `available_norms()`!");
  } // ... compute_norm(...)

  virtual std::vector< std::string > available_estimators() const override final
  {
    return {};
  }

  virtual double estimate(const VectorType& /*vector*/, const std::string /*type*/) const override final
  {
    std::abort();
  }

private:
  const size_t over_integrate_;
}; // class LinearEllipticEocStudy


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_EOCSTUDY_HH
