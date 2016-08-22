// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_INSTATIONARY_EOCSTUDY_HH
#define DUNE_GDT_TEST_INSTATIONARY_EOCSTUDY_HH

#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/provider/eoc.hh>

#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/discretizations/interfaces.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/prolongations.hh>

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \tparam ProblemType has to provide a type SolutionType which
 *         defines the type of the solution of the problem.
 * TODO: choose suitable SolutionType for Problems (provide Interface?)
 */
template <class GridImp, class ProblemImp>
class NonStationaryTestCase : public Stuff::Grid::Providers::EOC<GridImp>
{
  typedef Stuff::Grid::Providers::EOC<GridImp> EocBaseType;

public:
  typedef ProblemImp ProblemType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::SolutionType SolutionType;

public:
  template <class... Args>
  NonStationaryTestCase(const double divide_t_end_by_this, Args&&... args)
    : EocBaseType(std::forward<Args>(args)...)
    , divide_t_end_by_this_(divide_t_end_by_this)
    , zero_()
  {
  }

  virtual ~NonStationaryTestCase() = default;

  virtual const ProblemType& problem() const = 0;

  virtual void print_header(std::ostream& out = std::cout) const
  {
    out << "+===============================================================+\n"
        << "|+=============================================================+|\n"
        << "||  This is a GDT::Tests::NonStationaryTestCase, please provide ||\n"
        << "||  a meaningful message by implementing `print_header()`       ||\n"
        << "|+=============================================================+|\n"
        << "+===============================================================+" << std::endl;
  }

  virtual bool provides_exact_solution() const
  {
    return false;
  }

  virtual std::bitset<GridImp::dimension> periodic_directions() const
  {
    return std::bitset<GridImp::dimension>();
  }

  virtual double t_end() const
  {
    return problem().t_end() / divide_t_end_by_this_;
  }

  virtual const std::shared_ptr<const SolutionType> exact_solution() const
  {
    if (provides_exact_solution())
      DUNE_THROW(Stuff::Exceptions::you_have_to_implement_this,
                 "If provides_exact_solution() is true, exact_solution() has to be implemented!");
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call exact_solution() if provides_exact_solution() is false!");
    return zero_;
  }

private:
  const double divide_t_end_by_this_;
  const std::shared_ptr<const SolutionType> zero_;
}; // class NonStationaryTestCase


template <class TestCaseImp, class DiscretizerImp>
class NonStationaryEocStudy : public XT::Common::ConvergenceStudy
{
  typedef XT::Common::ConvergenceStudy BaseType;

protected:
  typedef TestCaseImp TestCaseType;
  typedef DiscretizerImp Discretizer;
  typedef typename Discretizer::DiscretizationType DiscretizationType;
  typedef typename TestCaseType::GridType GridType;
  typedef typename Discretizer::ProblemType ProblemType;
  typedef typename DiscretizationType::DiscreteSolutionType DiscreteSolutionType;
  typedef typename DiscretizationType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteSolutionType::key_type TimeFieldType;

  typedef typename TestCaseType::InitialValueType InitialValueType;
  typedef typename TestCaseType::template Level<Stuff::Grid::ChoosePartView::view>::Type GridViewType;

public:
  NonStationaryEocStudy(TestCaseType& test_case, const std::vector<std::string> only_these_norms = {},
                        const std::string visualize_prefix = "")
    : BaseType(only_these_norms)
    , test_case_(test_case)
    , current_refinement_(0)
    , last_computed_refinement_(std::numeric_limits<size_t>::max())
    , grid_widths_(num_refinements() + 1, -1.0)
    , time_to_solution_(0)
    , reference_solution_computed_(false)
    , discrete_exact_solution_computed_(false)
    , current_discretization_(nullptr)
    , current_solution_on_level_(nullptr)
    , reference_discretization_(nullptr)
    , reference_solution_(nullptr)
    , current_solution_(nullptr)
    , discrete_exact_solution_(nullptr)
    , visualize_prefix_(visualize_prefix)
  {
  }

  virtual ~NonStationaryEocStudy() = default;

  virtual size_t num_refinements() const override final
  {
    return test_case_.num_refinements();
  }

  virtual std::vector<std::string> provided_norms() const override final
  {
    std::vector<std::string> ret = available_norms();
    for (auto estimator : available_estimators()) {
      if (is_norm(estimator))
        DUNE_THROW(Stuff::Exceptions::internal_error,
                   "We do not want to handle the case that norms and estimators have the same name!");
      ret.push_back(estimator);
    }
    return ret;
  } // ... provided_norms(...)

  virtual double norm_reference_solution(const std::string type) override final
  {
    if (is_norm(type)) {
      compute_reference_solution();
      assert(reference_solution_);
      if (test_case_.provides_exact_solution()) {
        compute_discrete_exact_solution_vector();
        assert(discrete_exact_solution_);
        return compute_norm(*discrete_exact_solution_, type);
      } else {
        return compute_norm(*reference_solution_, type);
      }
    } else
      return 1.0;
  } // ... norm_reference_solution(...)

  virtual size_t current_num_DoFs() override final
  {
    assert(current_refinement_ <= num_refinements());
    const int level = test_case_.level_of(current_refinement_);
    return test_case_.grid().size(level, 0);
  } // ... current_num_DoFs(...)

  virtual size_t current_grid_size() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const int level = test_case_.level_of(current_refinement_);
    return test_case_.grid().size(level, 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() override final
  {
    assert(current_refinement_ <= num_refinements());
    if (grid_widths_[current_refinement_] < 0.0) {
      const int level      = test_case_.level_of(current_refinement_);
      const auto grid_view = test_case_.template level<Stuff::Grid::ChoosePartView::view>(level);
      Stuff::Grid::Dimensions<GridViewType> dimensions(grid_view);
      grid_widths_[current_refinement_] = dimensions.entity_width.max();
      assert(grid_widths_[current_refinement_] > 0.0);
    }
    return grid_widths_[current_refinement_];
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() override final
  {
    if (current_refinement_ != last_computed_refinement_) {
      assert(current_refinement_ <= num_refinements());
      // compute solution
      Timer timer;
      current_discretization_ = XT::Common::make_unique<DiscretizationType>(Discretizer::discretize(
          test_case_, test_case_, test_case_.level_of(current_refinement_), test_case_.periodic_directions()));
      current_solution_on_level_ = XT::Common::make_unique<DiscreteSolutionType>(current_discretization_->solve());
      time_to_solution_          = timer.elapsed();
      // prolong to reference grid part
      compute_reference_solution();
      assert(reference_solution_);
      if (!current_solution_)
        current_solution_ = XT::Common::make_unique<DiscreteSolutionType>(*reference_solution_);
      // time prolongation
      DiscreteSolutionType time_prolongated_current_solution_on_level;
      const auto time_prolongated_current_solution_on_level_it_end = time_prolongated_current_solution_on_level.end();
      auto current_solution_on_level_it                            = current_solution_on_level_->begin();
      const auto current_solution_on_level_it_last                 = --current_solution_on_level_->end();
      const auto current_solution_it_end                           = current_solution_->end();
      auto last_time                                               = current_solution_->begin()->first;
      for (auto current_solution_it = current_solution_->begin(); current_solution_it != current_solution_it_end;
           ++current_solution_it) {
        const auto time          = current_solution_it->first;
        const auto time_on_level = current_solution_on_level_it->first;
        const auto inserted_it   = time_prolongated_current_solution_on_level.emplace_hint(
            time_prolongated_current_solution_on_level_it_end, time, current_solution_on_level_it->second);
        if (time_on_level < time && current_solution_on_level_it != current_solution_on_level_it_last) {
          // compute weighted average of the two values of current_solution_on_level_
          inserted_it->second.vector() = (current_solution_on_level_it->second.vector() * (time_on_level - last_time)
                                          + (++current_solution_on_level_it)->second.vector() * (time - time_on_level))
                                         * (1.0 / (time - last_time));
        }
        last_time = time;
      }
      // spatial prolongation
      for (auto current_solution_it = current_solution_->begin(); current_solution_it != current_solution_it_end;
           ++current_solution_it)
        prolong(time_prolongated_current_solution_on_level.at(current_solution_it->first),
                (*current_solution_it).second);
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        size_t counter = 0;
        for (auto current_solution_it = current_solution_->begin(); current_solution_it != current_solution_it_end;
             ++current_solution_it, ++counter)
          current_solution_it->second.visualize(visualize_prefix_ + "_solution_"
                                                    + Dune::XT::Common::to_string(current_refinement_),
                                                Dune::XT::Common::to_string(counter));
      }
    }
    return time_to_solution_;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override final
  {
    // get current solution
    assert(current_refinement_ <= num_refinements());
    compute_on_current_refinement();
    assert(last_computed_refinement_ == current_refinement_);
    if (is_norm(type)) {
      assert(current_solution_);
      compute_reference_solution();
      assert(reference_discretization_);
      // compute error
      std::unique_ptr<DiscreteSolutionType> difference;
      if (test_case_.provides_exact_solution()) {
        compute_discrete_exact_solution_vector();
        assert(discrete_exact_solution_);
        difference = Dune::XT::Common::make_unique<DiscreteSolutionType>(*discrete_exact_solution_);
      } else {
        assert(reference_solution_);
        difference = Dune::XT::Common::make_unique<DiscreteSolutionType>(*reference_solution_);
      }
      for (auto& difference_at_time : *difference) {
        auto time = difference_at_time.first;
        assert(current_solution_->count(time) && "Time steps must be the same");
        difference_at_time.second.vector() = difference_at_time.second.vector() - current_solution_->at(time).vector();
      }
      return compute_norm(*difference, type);
    } else {
      assert(current_solution_on_level_);
      return estimate(*current_solution_on_level_, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override final
  {
    if (current_refinement_ <= num_refinements())
      ++current_refinement_;
  }

  std::map<std::string, std::vector<double>> run(std::ostream& out, const bool print_timings = false)
  {
    return BaseType::run(false, out, print_timings);
  }

protected:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_ = XT::Common::make_unique<DiscretizationType>(Discretizer::discretize(
          test_case_, test_case_, test_case_.reference_level(), test_case_.periodic_directions()));
      reference_solution_ = XT::Common::make_unique<DiscreteSolutionType>(reference_discretization_->solve());
      reference_solution_computed_ = true;
      if (!visualize_prefix_.empty()) {
        size_t counter                       = 0;
        const auto reference_solution_it_end = reference_solution_->end();
        for (auto reference_solution_it = reference_solution_->begin();
             reference_solution_it != reference_solution_it_end;
             ++reference_solution_it, ++counter)
          reference_solution_it->second.visualize(visualize_prefix_ + "_reference"
                                                      + Dune::XT::Common::to_string(current_refinement_),
                                                  Dune::XT::Common::to_string(counter));
      }
    }
  } // ... compute_reference_solution()

  void compute_discrete_exact_solution_vector()
  {
    if (!discrete_exact_solution_computed_) {
      discrete_exact_solution_ = Dune::XT::Common::make_unique<DiscreteSolutionType>();
      compute_reference_solution();
      const auto exact_solution            = test_case_.exact_solution();
      const auto reference_solution_it_end = reference_solution_->end();
      for (auto reference_solution_it = reference_solution_->begin();
           reference_solution_it != reference_solution_it_end;
           ++reference_solution_it) {
        const double time                          = reference_solution_it->first;
        const auto discrete_exact_solution_at_time = exact_solution->evaluate_at_time(time);
        const auto inserted_it                     = discrete_exact_solution_->emplace_hint(
            discrete_exact_solution_->end(), time, reference_solution_it->second);
        project(*discrete_exact_solution_at_time, inserted_it->second);
      }
      if (!visualize_prefix_.empty()) {
        size_t counter                            = 0;
        const auto discrete_exact_solution_it_end = discrete_exact_solution_->end();
        for (auto discrete_exact_solution_it = discrete_exact_solution_->begin();
             discrete_exact_solution_it != discrete_exact_solution_it_end;
             ++discrete_exact_solution_it, ++counter)
          discrete_exact_solution_it->second.visualize(visualize_prefix_ + "_exact_solution"
                                                           + Dune::XT::Common::to_string(current_refinement_),
                                                       Dune::XT::Common::to_string(counter));
      }
      discrete_exact_solution_computed_ = true;
    }
  }

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  }

  virtual std::vector<std::string> available_norms() const = 0;

  virtual std::vector<std::string> available_estimators() const = 0;

  virtual double estimate(const DiscreteSolutionType& solution, const std::string type) = 0;

  virtual double compute_norm(const DiscreteSolutionType& solution, const std::string type) = 0;

  TestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  mutable std::vector<double> grid_widths_;
  double time_to_solution_;
  bool reference_solution_computed_;
  bool discrete_exact_solution_computed_;
  std::unique_ptr<DiscretizationType> current_discretization_;
  std::unique_ptr<DiscreteSolutionType> current_solution_on_level_;
  std::unique_ptr<DiscretizationType> reference_discretization_;
  std::unique_ptr<DiscreteSolutionType> reference_solution_;
  std::unique_ptr<DiscreteSolutionType> current_solution_;
  std::unique_ptr<DiscreteSolutionType> discrete_exact_solution_;
  const std::string visualize_prefix_;
}; // class NonStationaryEocStudy


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDY_HH
