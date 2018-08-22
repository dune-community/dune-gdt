// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TEST_STATIONARY_EOCSTUDY_HH
#define DUNE_GDT_TEST_STATIONARY_EOCSTUDY_HH

#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/information.hh>

#include <dune/gdt/prolongations.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class TestCaseImp, class DiscretizerImp>
class StationaryEocStudy : public XT::Common::ConvergenceStudy
{
  typedef XT::Common::ConvergenceStudy BaseType;

protected:
  typedef TestCaseImp TestCaseType;
  typedef DiscretizerImp Discretizer;
  typedef typename Discretizer::DiscretizationType DiscretizationType;
  typedef typename TestCaseType::GridType GridType;
  typedef typename Discretizer::ProblemType ProblemType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::AnsatzSpaceType SpaceType;
  typedef GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;

  typedef typename TestCaseType::FunctionType FunctionType;
  typedef typename TestCaseType::LevelGridViewType GridViewType;

public:
  StationaryEocStudy(TestCaseType& test_case,
                     const std::vector<std::string> only_these_norms = {},
                     const std::string visualize_prefix = "",
                     XT::Common::Configuration solver_options = DiscretizationType::solver_options())
    : BaseType(only_these_norms)
    , test_case_(test_case)
    , current_refinement_(0)
    , last_computed_refinement_(std::numeric_limits<size_t>::max())
    , grid_widths_(test_case.num_refinements() + 1, -1.0)
    , time_to_solution_(0)
    , reference_solution_computed_(false)
    , current_discretization_(nullptr)
    , current_solution_vector_on_level_(nullptr)
    , reference_discretization_(nullptr)
    , reference_solution_vector_(nullptr)
    , current_solution_vector_(nullptr)
    , visualize_prefix_(visualize_prefix)
    , current_num_DoFs_(0)
    , solver_options_(solver_options)
  {
  }

  virtual ~StationaryEocStudy() = default;

  virtual size_t num_refinements() const override final
  {
    return test_case_.num_refinements();
  }

  virtual std::vector<std::string> provided_norms() const override final
  {
    std::vector<std::string> ret = available_norms();
    for (auto estimator : available_estimators()) {
      if (is_norm(estimator))
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "We do not want to handle the case that norms and estimators have the same name!");
      ret.push_back(estimator);
    }
    return ret;
  } // ... provided_norms(...)

  virtual double norm_reference_solution(const std::string type) override final
  {
    if (!is_norm(type))
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Do not call norm_reference_solution() for an estimator!\n"
                     << "type: "
                     << type
                     << "\n");
    if (test_case_.provides_exact_solution()) {
      // visualize
      if (!visualize_prefix_.empty()) {
        test_case_.exact_solution().visualize(test_case_.reference_grid_view(), visualize_prefix_ + "_exact_solution");
      }
      return compute_norm(test_case_.reference_grid_view(), test_case_.exact_solution(), type);
    } else {
      compute_reference_solution();
      DXT_ASSERT(reference_discretization_);
      DXT_ASSERT(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->ansatz_space(), *reference_solution_vector_, "reference solution");
      return compute_norm(test_case_.reference_grid_view(), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_num_DoFs() override final
  {
    if (current_refinement_ != last_computed_refinement_) {
      DXT_ASSERT(current_refinement_ <= num_refinements());
      current_num_DoFs_ = Discretizer::discretize(test_case_.level_provider(current_refinement_),
                                                  test_case_.problem(),
                                                  test_case_.level_of(current_refinement_))
                              .ansatz_space()
                              .mapper()
                              .size();
    }
    return current_num_DoFs_;
  } // ... current_num_DoFs(...)

  virtual size_t current_grid_size() const override final
  {
    DXT_ASSERT(current_refinement_ <= num_refinements());
    return test_case_.level_provider(current_refinement_)
        .template layer<TestCaseType::layer_type, XT::Grid::Backends::view>(test_case_.level_of(current_refinement_))
        .indexSet()
        .size(0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() override final
  {
    DXT_ASSERT(current_refinement_ <= num_refinements());
    if (grid_widths_[current_refinement_] < 0.0) {
      const int level = test_case_.level_of(current_refinement_);
      const auto grid_layer = test_case_.level_provider(current_refinement_)
                                  .template layer<TestCaseType::layer_type, XT::Grid::Backends::view>(level);
      grid_widths_[current_refinement_] = XT::Grid::dimensions(grid_layer).entity_width.max();
      DXT_ASSERT(grid_widths_[current_refinement_] > 0.0);
    }
    return grid_widths_[current_refinement_];
  } // ... current_grid_width(...)


  virtual double compute_on_current_refinement() override final
  {
    if (current_refinement_ != last_computed_refinement_) {
      DXT_ASSERT(current_refinement_ <= num_refinements());
      // compute solution
      Timer timer;
      current_discretization_ = XT::Common::make_unique<DiscretizationType>(
          Discretizer::discretize(test_case_.level_provider(current_refinement_),
                                  test_case_.problem(),
                                  test_case_.level_of(current_refinement_)));
      current_solution_vector_on_level_ =
          XT::Common::make_unique<VectorType>(current_discretization_->solve(solver_options_));
      time_to_solution_ = timer.elapsed();
      const ConstDiscreteFunctionType current_refinement_solution(
          current_discretization_->ansatz_space(), *current_solution_vector_on_level_, "solution on current level");
      // prolong to reference grid part
      compute_reference_solution();
      DXT_ASSERT(reference_discretization_);
      if (!current_solution_vector_)
        current_solution_vector_ = XT::Common::make_unique<VectorType>(reference_discretization_->create_vector());
      DiscreteFunctionType reference_refinement_solution(
          reference_discretization_->ansatz_space(), *current_solution_vector_, "solution on reference grid part");
      prolong(current_refinement_solution, reference_refinement_solution);
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(test_case_.reference_grid_view(),
                                             visualize_prefix_ + "_problem_"
                                                 + Dune::XT::Common::to_string(current_refinement_));
        current_refinement_solution.visualize(visualize_prefix_ + "_solution_"
                                              + Dune::XT::Common::to_string(current_refinement_));
      }
    }
    return time_to_solution_;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override final
  {
    // get current solution
    DXT_ASSERT(current_refinement_ <= num_refinements());
    compute_on_current_refinement();
    DXT_ASSERT(last_computed_refinement_ == current_refinement_);
    if (is_norm(type)) {
      DXT_ASSERT(current_solution_vector_);
      compute_reference_solution();
      DXT_ASSERT(reference_discretization_);
      const ConstDiscreteFunctionType current_solution(
          reference_discretization_->ansatz_space(), *current_solution_vector_, "current solution");
      // compute error
      if (test_case_.provides_exact_solution()) {
        return compute_norm(test_case_.reference_grid_view(), test_case_.exact_solution() - current_solution, type);
      } else {
        // get reference solution
        compute_reference_solution();
        DXT_ASSERT(reference_discretization_);
        DXT_ASSERT(reference_solution_vector_);
        const ConstDiscreteFunctionType reference_solution(
            reference_discretization_->ansatz_space(), *reference_solution_vector_, "reference solution");
        return compute_norm(test_case_.reference_grid_view(), reference_solution - current_solution, type);
      }
    } else {
      DXT_ASSERT(current_solution_vector_on_level_);
      return estimate(*current_solution_vector_on_level_, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override final
  {
    if (current_refinement_ <= num_refinements())
      ++current_refinement_;
  }

  std::map<std::string, std::vector<double>> run(std::ostream& out, const bool print_timings = true)
  {
    return BaseType::run(true, out, print_timings);
  }

protected:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_ = XT::Common::make_unique<DiscretizationType>(
          Discretizer::discretize(test_case_.reference_provider(), test_case_.problem(), test_case_.reference_level()));
      reference_solution_vector_ =
          XT::Common::make_unique<VectorType>(reference_discretization_->solve(solver_options_));

      reference_solution_computed_ = true;
      // visualize
      if (!visualize_prefix_.empty()) {
        this->test_case_.problem().visualize(test_case_.reference_grid_view(),
                                             visualize_prefix_ + "_problem_reference");
        ConstDiscreteFunctionType(
            reference_discretization_->ansatz_space(), *reference_solution_vector_, "reference solution")
            .visualize(visualize_prefix_ + "_reference_solution");
      }
    }
  } // ... compute_reference_solution()

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  }

  virtual std::vector<std::string> available_norms() const = 0;

  virtual std::vector<std::string> available_estimators() const = 0;

  virtual double estimate(const VectorType& vector, const std::string type) = 0;

  virtual double compute_norm(const GridViewType& grid_view, const FunctionType& function, const std::string type) = 0;

  TestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  std::vector<double> grid_widths_;
  double time_to_solution_;
  bool reference_solution_computed_;
  std::unique_ptr<DiscretizationType> current_discretization_;
  std::unique_ptr<VectorType> current_solution_vector_on_level_;
  std::unique_ptr<DiscretizationType> reference_discretization_;
  std::unique_ptr<VectorType> reference_solution_vector_;
  std::unique_ptr<VectorType> current_solution_vector_;
  const std::string visualize_prefix_;
  size_t current_num_DoFs_;
  const XT::Common::Configuration solver_options_;
}; // class StationaryEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_EOCSTUDY_HH
