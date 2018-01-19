// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATORS_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATORS_HH

#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/common/test/common.hh>

#include "discretizers/ipdg.hh"
#include "eocstudy.hh"
#include "estimators/swipdg-fluxreconstruction.hh"
#include "swipdg-estimator-expectations.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class TestCaseImp, class DiscretizerImp>
class LinearEllipticSwipdgEstimatorStudy : public LinearEllipticEocStudy<TestCaseImp, DiscretizerImp>
{
  typedef LinearEllipticSwipdgEstimatorStudy<TestCaseImp, DiscretizerImp> ThisType;
  typedef LinearEllipticEocStudy<TestCaseImp, DiscretizerImp> BaseType;

public:
  using typename BaseType::TestCaseType;
  using typename BaseType::Discretizer;
  using typename BaseType::DiscretizationType;
  using typename BaseType::FunctionType;
  using typename BaseType::VectorType;
  using typename BaseType::SpaceType;
  using typename BaseType::ProblemType;
  using typename BaseType::GridViewType;

private:
  typedef typename ProblemType::DiffusionFactorType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType DiffusionTensorType;

  static const int polOrder = Discretizer::polOrder;

  typedef LinearElliptic::SwipdgFluxreconstrutionEstimators::
      LocalNonconformityESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
          LocalNonconformityESV2007Estimator;
  typedef LinearElliptic::SwipdgFluxreconstrutionEstimators::
      LocalResidualESV2007<SpaceType, VectorType, FunctionType, DiffusionFactorType, DiffusionTensorType, GridViewType>
          LocalResidualESV2007Estimator;
  typedef LinearElliptic::SwipdgFluxreconstrutionEstimators::
      LocalDiffusiveFluxESV2007<SpaceType, VectorType, DiffusionFactorType, DiffusionTensorType, GridViewType>
          LocalDiffusiveFluxESV2007Estimator;
  typedef LinearElliptic::SwipdgFluxreconstrutionEstimators::
      ESV2007<SpaceType, VectorType, FunctionType, DiffusionFactorType, DiffusionTensorType, GridViewType>
          ESV2007Estimator;
  typedef LinearElliptic::SwipdgFluxreconstrutionEstimators::ESV2007AlternativeSummation<SpaceType,
                                                                                         VectorType,
                                                                                         FunctionType,
                                                                                         DiffusionFactorType,
                                                                                         DiffusionTensorType,
                                                                                         GridViewType>
      ESV2007AlternativeSummationEstimator;

public:
  // a perfect forwarding ctor did not do the job here, since it was not able to match the std::initializer_list: {"L2"}
  LinearEllipticSwipdgEstimatorStudy(TestCaseType& test_case,
                                     const std::vector<std::string> only_these_norms = {},
                                     const std::string visualize_prefix = "",
                                     const size_t over_integrate = 2)
    : BaseType(test_case, only_these_norms, visualize_prefix, over_integrate)
  {
  }

  virtual ~LinearEllipticSwipdgEstimatorStudy() = default;

  virtual std::string identifier() const override final
  {
    return "gdt.linearelliptic.estimators.swipdg.polorder_" + Dune::XT::Common::to_string(int(polOrder));
  }

  virtual size_t expected_rate(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker, see the explanation in LinearEllipticEocStudy!
    return LinearEllipticSwipdgEstimatorExpectations<TestCaseType, Discretizer::type, polOrder>::rate(type);
  }

  virtual std::vector<double> expected_results(const std::string type) const override final
  {
    // If you get an undefined reference here from the linker, see above!
    return LinearEllipticSwipdgEstimatorExpectations<TestCaseType, Discretizer::type, polOrder>::results(
        this->test_case_, type);
  }

  virtual std::vector<std::string> available_norms() const override final
  {
    return {"energy"};
  }

  virtual std::vector<std::string> available_estimators() const override final
  {
    return {LocalNonconformityESV2007Estimator::id(),
            LocalResidualESV2007Estimator::id(),
            LocalDiffusiveFluxESV2007Estimator::id(),
            "efficiency_" + ESV2007Estimator::id(),
            "efficiency_" + ESV2007AlternativeSummationEstimator::id()};
  }

  virtual double estimate(const VectorType& vector, const std::string type) override final
  {
    const auto& space = this->current_discretization_->ansatz_space();
    const auto& force = this->test_case_.problem().force();
    const auto& diffusion_factor = this->test_case_.problem().diffusion_factor();
    const auto& diffusion_tensor = this->test_case_.problem().diffusion_tensor();
    const auto grid_view = this->test_case_.level_provider(this->current_refinement_)
                               .template layer<TestCaseType::layer_type, XT::Grid::Backends::view>(
                                   this->test_case_.level_of(this->current_refinement_));
    if (type == LocalNonconformityESV2007Estimator::id())
      return LocalNonconformityESV2007Estimator::estimate(grid_view, space, vector, diffusion_factor, diffusion_tensor);
    else if (type == LocalResidualESV2007Estimator::id())
      return LocalResidualESV2007Estimator::estimate(
          grid_view, space, vector, force, diffusion_factor, diffusion_tensor);
    else if (type == LocalDiffusiveFluxESV2007Estimator::id())
      return LocalDiffusiveFluxESV2007Estimator::estimate(grid_view, space, vector, diffusion_factor, diffusion_tensor);
    else if (type == ESV2007Estimator::id())
      return ESV2007Estimator::estimate(grid_view, space, vector, force, diffusion_factor, diffusion_tensor);
    else if (type == "efficiency_" + ESV2007Estimator::id())
      return estimate(vector, ESV2007Estimator::id()) / this->current_error_norm("energy");
    else if (type == ESV2007AlternativeSummationEstimator::id())
      return ESV2007AlternativeSummationEstimator::estimate(
          grid_view, space, vector, force, diffusion_factor, diffusion_tensor);
    else if (type == "efficiency_" + ESV2007AlternativeSummationEstimator::id())
      return estimate(vector, ESV2007AlternativeSummationEstimator::id()) / this->current_error_norm("energy");
    else
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                 "Wrong type `" << type << "` requested (see available_estimators()!");
    return 0.;
  } // ... estimate(...)

  std::map<std::string, std::vector<double>> run(std::ostream& out, const bool print_timings = false)
  {
    return XT::Common::ConvergenceStudy::run(false, out, print_timings);
  }
}; // class LinearEllipticSwipdgEstimatorStudy


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class TestCaseType>
struct linearelliptic_SWIPDG_estimators : public ::testing::Test
{
  template <Dune::GDT::Backends space_backend, Dune::XT::LA::Backends la_backend, int polOrder>
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
    TestCaseType test_case;
    test_case.print_header(DXTC_LOG_INFO_0);
    DXTC_LOG_INFO_0 << std::endl;
    typedef LinearElliptic::IpdgDiscretizer<typename TestCaseType::GridType,
                                            TestCaseType::layer_type,
                                            space_backend,
                                            la_backend,
                                            polOrder,
                                            typename TestCaseType::ProblemType::RangeFieldType,
                                            1,
                                            LocalEllipticIpdgIntegrands::Method::swipdg>
        Discretizer;
    Dune::GDT::Test::LinearEllipticSwipdgEstimatorStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(
          eoc_study, eoc_study.run(DXTC_LOG_INFO_0), /*zero_tolerance=*/1.4e-10);
    } catch (Dune::XT::Common::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.swipdg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_SWIPDG_estimators


#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_SWIPDG_ESTIMATORS_HH
