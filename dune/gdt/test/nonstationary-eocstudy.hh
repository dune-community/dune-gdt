// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_NONSTATIONARY_EOCSTUDY_HH
#define DUNE_GDT_TEST_NONSTATIONARY_EOCSTUDY_HH

#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/provider/eoc.hh>

#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/discretizations/interfaces.hh>
#include <dune/gdt/prolongations.hh>

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \tparam ProblemType has to provide a type SolutionType which
 *         defines the type of the solution of the problem.
 * TODO: choose suitable SolutionType for Problems (provide Interface?)
 */
template< class GridImp, class ProblemImp >
class NonStationaryTestCase
  : public Stuff::Grid::Providers::EOC< GridImp >
{
  typedef Stuff::Grid::Providers::EOC< GridImp > EocBaseType;
public:
  typedef ProblemImp                         ProblemType;
  typedef typename ProblemType::InitialValueType InitialValueType;
  typedef typename ProblemType::SolutionType SolutionType;

public:
  template< class... Args >
  NonStationaryTestCase(Args&& ...args)
    : EocBaseType(std::forward< Args >(args)...)
    , zero_()
  {}

  virtual ~NonStationaryTestCase() {}

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

  virtual std::bitset< GridImp::dimension > periodic_directions() const
  {
    return std::bitset< GridImp::dimension >();
  }

  virtual const std::shared_ptr< const SolutionType > exact_solution() const
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
  const std::shared_ptr< const SolutionType > zero_;
}; // class NonStationaryTestCase


template< class TestCaseImp, class DiscretizerImp >
class NonStationaryEocStudy
  : public Stuff::Common::ConvergenceStudy
{
  typedef Stuff::Common::ConvergenceStudy          BaseType;
protected:
  typedef TestCaseImp                              TestCaseType;
  typedef DiscretizerImp                           Discretizer;
  typedef typename Discretizer::DiscretizationType DiscretizationType;
  typedef typename TestCaseType::GridType          GridType;
  typedef typename Discretizer::ProblemType        ProblemType;
  typedef typename DiscretizationType::DiscreteSolutionType  DiscreteSolutionType;
  typedef typename DiscretizationType::DiscreteFunctionType DiscreteFunctionType;

  typedef typename TestCaseType::InitialValueType InitialValueType;
  typedef typename TestCaseType::template Level< Stuff::Grid::ChoosePartView::view >::Type GridViewType;

public:
  NonStationaryEocStudy(TestCaseType& test_case,
                     const std::vector< std::string > only_these_norms = {},
                     const std::string visualize_prefix = "")
    : BaseType(only_these_norms)
    , test_case_(test_case)
    , current_refinement_(0)
    , last_computed_refinement_(std::numeric_limits< size_t >::max())
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
  {}

  virtual ~NonStationaryEocStudy() {}

  virtual size_t num_refinements() const override final
  {
    return test_case_.num_refinements();
  }

  virtual std::vector< std::string > provided_norms() const override final
  {
    std::vector< std::string > ret = available_norms();
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

  virtual size_t current_grid_size() const override final
  {
    assert(current_refinement_ <= num_refinements());
    const int level = test_case_.level_of(current_refinement_);
    return test_case_.grid().size(level, 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const override final
  {
    assert(current_refinement_ <= num_refinements());
    if (grid_widths_[current_refinement_] < 0.0) {
      const int level = test_case_.level_of(current_refinement_);
      const auto grid_view = test_case_.template level< Stuff::Grid::ChoosePartView::view >(level);
      Stuff::Grid::Dimensions< GridViewType > dimensions(grid_view);
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
      current_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(Discretizer::discretize(test_case_,
                                                                                     test_case_.problem(),
                                                                                     test_case_.level_of(current_refinement_),
                                                                                     test_case_.periodic_directions()));
      current_solution_on_level_ = Stuff::Common::make_unique< DiscreteSolutionType >(current_discretization_->solve(test_case_.problem().is_linear()));
      time_to_solution_ = timer.elapsed();
      // prolong to reference grid part
      compute_reference_solution();
      assert(reference_solution_);
      if (!current_solution_)
        current_solution_ = Stuff::Common::make_unique< DiscreteSolutionType >(*reference_solution_);
      // time prolongation
      DiscreteSolutionType time_prolongated_current_solution_on_level;
      size_t jj = 0;
      for (size_t ii = 0; ii < current_solution_->size(); ++ii) {
        const auto time_ii = current_solution_->operator[](ii).first;
        const auto time_jj = current_solution_on_level_->operator[](jj).first;
        time_prolongated_current_solution_on_level.emplace_back(std::make_pair(time_ii, current_solution_on_level_->operator[](jj).second));
        if (time_jj < time_ii && jj != current_solution_on_level_->size() - 1) {
          // compute weighted average of the two values (at time jj and jj+1) of current_solution_on_level_
          const auto time_ii_minus_1 = current_solution_->operator[](ii-1).first;
          time_prolongated_current_solution_on_level[ii].second.vector() = (current_solution_on_level_->operator[](jj).second.vector()*(time_jj - time_ii_minus_1)
                                          + current_solution_on_level_->operator[](jj+1).second.vector()*(time_ii - time_jj))
                                         *(1.0/(time_ii - time_ii_minus_1));
          ++jj;
        }
      }
      // spatial prolongation
      for (size_t ii = 0; ii < current_solution_->size(); ++ii)
        Operators::prolong(time_prolongated_current_solution_on_level[ii].second, current_solution_->operator[](ii).second);
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        for (size_t ii = 0; ii < current_solution_->size(); ++ii) {
          current_solution_->operator[](ii).second.visualize(visualize_prefix_ + "_solution_" + DSC::to_string(current_refinement_), DSC::to_string(ii));
        }
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
      if (test_case_.provides_exact_solution()) {
        compute_discrete_exact_solution_vector();
        assert(discrete_exact_solution_);
        std::unique_ptr< DiscreteSolutionType > difference = DSC::make_unique< DiscreteSolutionType >(*discrete_exact_solution_);
        for (size_t ii = 0; ii < difference->size(); ++ii) {
          assert(DSC::FloatCmp::eq(difference->operator[](ii).first, current_solution_->operator[](ii).first) && "Time steps must be the same");
          difference->operator[](ii).second.vector() = difference->operator[](ii).second.vector() - current_solution_->operator[](ii).second.vector();
        }
        return compute_norm(*difference, type);
      } else {
        assert(reference_solution_);
        std::unique_ptr< DiscreteSolutionType > difference = DSC::make_unique< DiscreteSolutionType >(*reference_solution_);
        for (size_t ii = 0; ii < difference->size(); ++ii) {
          assert(DSC::FloatCmp::eq(difference->operator[](ii).first, current_solution_->operator[](ii).first) && "Time steps must be the same");
          difference->operator[](ii).second.vector() = difference->operator[](ii).second.vector() - current_solution_->operator[](ii).second.vector();
        }
        return compute_norm(*difference, type);
      }
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

  std::map< std::string, std::vector< double > > run(std::ostream& out, const bool print_timings = false)
  {
    return BaseType::run(false, out, print_timings);
  }

protected:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_
          = Stuff::Common::make_unique< DiscretizationType >(Discretizer::discretize(test_case_,
                                                                                     test_case_.problem(),
                                                                                     test_case_.reference_level(),
                                                                                     test_case_.periodic_directions()));
      reference_solution_ = Stuff::Common::make_unique< DiscreteSolutionType >(reference_discretization_->solve(test_case_.problem().is_linear()));
      reference_solution_computed_ = true;
      if (!visualize_prefix_.empty()) {
        for (size_t ii = 0; ii < reference_solution_->size(); ++ii) {
          reference_solution_->operator[](ii).second.visualize(visualize_prefix_ + "_reference", DSC::to_string(ii));
        }
      }
    }
  } // ... compute_reference_solution()

  void compute_discrete_exact_solution_vector()
  {
    if (!discrete_exact_solution_computed_) {
      discrete_exact_solution_ = DSC::make_unique< DiscreteSolutionType >();
      compute_reference_solution();
      const auto exact_solution = test_case_.exact_solution();
      for (size_t ii = 0; ii < reference_solution_->size(); ++ii) {
        const double time = reference_solution_->operator[](ii).first;
        const auto discrete_exact_solution_at_time = exact_solution->evaluate_at_time(time);
        (*discrete_exact_solution_).emplace_back(std::make_pair(time, DiscreteFunctionType(reference_discretization_->fv_space(), "exact solution at time " + DSC::to_string(time))));
        project(*discrete_exact_solution_at_time, discrete_exact_solution_->operator[](ii).second);
      }
      if (!visualize_prefix_.empty()) {
        for (size_t ii = 0; ii < discrete_exact_solution_->size(); ++ii) {
          discrete_exact_solution_->operator[](ii).second.visualize(visualize_prefix_ + "_exact_solution", DSC::to_string(ii));
        }
      }
      discrete_exact_solution_computed_ = true;
    }
  }

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  }

  virtual std::vector< std::string > available_norms() const = 0;

  virtual std::vector< std::string > available_estimators() const = 0;

  virtual double estimate(const DiscreteSolutionType& solution, const std::string type) const = 0;

  virtual double compute_norm(const DiscreteSolutionType& solution,
                              const std::string type) const = 0;

  TestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  mutable std::vector< double > grid_widths_;
  double time_to_solution_;
  bool reference_solution_computed_;
  bool discrete_exact_solution_computed_;
  std::unique_ptr< DiscretizationType > current_discretization_;
  std::unique_ptr< DiscreteSolutionType > current_solution_on_level_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
  std::unique_ptr< DiscreteSolutionType > reference_solution_;
  std::unique_ptr< DiscreteSolutionType > current_solution_;
  std::unique_ptr< DiscreteSolutionType > discrete_exact_solution_;
  const std::string visualize_prefix_;
}; // class NonStationaryEocStudy


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_NONSTATIONARY_EOCSTUDY_HH
