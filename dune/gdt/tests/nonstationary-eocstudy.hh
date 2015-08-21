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
#include <dune/gdt/operators/prolongations.hh>

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \tparam ProblemType has to provide a type FunctionType (derived from Stuff::LocalizableFunctionInterface) which
 *         defines the type of the solution of the problem.
 */
template< class GridImp, class ProblemImp >
class NonStationaryTestCase
  : public Stuff::Grid::Providers::EOC< GridImp >
{
  typedef Stuff::Grid::Providers::EOC< GridImp > EocBaseType;
public:
  typedef ProblemImp                         ProblemType;
  typedef typename ProblemType::FunctionType FunctionType;
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
  typedef typename DiscretizationType::VectorType  VectorType;
  typedef typename DiscretizationType::StationaryVectorType StationaryVectorType;
  typedef typename DiscretizationType::FVSpaceType SpaceType;
  typedef GDT::DiscreteFunction< SpaceType, StationaryVectorType >      DiscreteFunctionType;
  typedef GDT::ConstDiscreteFunction< SpaceType, StationaryVectorType > ConstDiscreteFunctionType;

  typedef typename TestCaseType::FunctionType FunctionType;
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
    , exact_solution_vector_computed_(false)
    , current_discretization_(nullptr)
    , current_solution_vector_on_level_(nullptr)
    , reference_discretization_(nullptr)
    , reference_solution_vector_(nullptr)
    , current_solution_vector_(nullptr)
    , exact_solution_vector_(nullptr)
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
      assert(reference_solution_vector_);
      if (test_case_.provides_exact_solution()) {
        compute_exact_solution_vector();
        assert(exact_solution_vector_);
        return compute_norm(test_case_.reference_grid_view(), *exact_solution_vector_, type);
      } else {
        return compute_norm(test_case_.reference_grid_view(), *reference_solution_vector_, type);
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
                                                                                     test_case_.level_of(current_refinement_)));
      current_solution_vector_on_level_ = Stuff::Common::make_unique< VectorType >(current_discretization_->solve(test_case_.problem().is_linear()));
      time_to_solution_ = timer.elapsed();
      // prolong to reference grid part
      compute_reference_solution();
      assert(reference_solution_vector_);
      if (!current_solution_vector_)
        current_solution_vector_ = Stuff::Common::make_unique< VectorType >(*reference_solution_vector_);
      size_t jj = 0;
      const auto current_solution_vector_on_level_copy = DSC::make_unique< VectorType >(*current_solution_vector_on_level_);
      current_solution_vector_on_level_->clear();
      // time prolongation
      for (size_t ii = 0; ii < current_solution_vector_->size(); ++ii) {
        const auto& curr_sol_on_level_jj = current_solution_vector_on_level_copy->operator[](jj);
        const auto& curr_sol_ii = current_solution_vector_->operator[](ii);
        if (DSC::FloatCmp::ge(curr_sol_on_level_jj.first, curr_sol_ii.first) || jj == current_solution_vector_on_level_copy->size() - 1) {
          current_solution_vector_on_level_->emplace_back(std::make_pair(curr_sol_ii.first, curr_sol_on_level_jj.second));
        } else {
          // compute average solution
          const auto& curr_sol_on_level_jj_plus_1 = current_solution_vector_on_level_copy->operator[](jj+1);
          typename VectorType::value_type average_pair = curr_sol_on_level_jj;
          average_pair.first = curr_sol_ii.first;
          const auto& curr_sol_ii_minus_1 = current_solution_vector_->operator[](ii-1);
          average_pair.second = (curr_sol_on_level_jj.second*(curr_sol_on_level_jj.first - curr_sol_ii_minus_1.first)
                                          + curr_sol_on_level_jj_plus_1.second*(curr_sol_ii.first - curr_sol_on_level_jj.first))
                                         *(1.0/(curr_sol_ii.first - curr_sol_ii_minus_1.first));
          current_solution_vector_on_level_->emplace_back(average_pair);
          ++jj;
        }
      }
      // spatial prolongation
      for (size_t ii = 0; ii < current_solution_vector_->size(); ++ii) {
        DiscreteFunctionType curr_ii_as_discrete_func(current_discretization_->fv_space(),
                                                      current_solution_vector_on_level_->operator[](ii).second,
                                                      "solution on current level");
        DiscreteFunctionType curr_ii_on_ref_as_discrete_func(reference_discretization_->fv_space(),
                                                             current_solution_vector_->operator[](ii).second,
                                                             "solution on reference level");
        Operators::prolong(curr_ii_as_discrete_func, curr_ii_on_ref_as_discrete_func);
      }
      last_computed_refinement_ = current_refinement_;
      // visualize
      if (!visualize_prefix_.empty()) {
        DiscreteFunctionType curr_ii_on_ref_as_discrete_func(reference_discretization_->fv_space(),
                                                             current_solution_vector_->operator[](0).second,
                                                             "solution on reference level");
        for (size_t ii = 0; ii < current_solution_vector_->size(); ++ii) {
          curr_ii_on_ref_as_discrete_func.vector() = current_solution_vector_->operator[](ii).second;
          curr_ii_on_ref_as_discrete_func.template visualize_factor< 0 >(visualize_prefix_ + "_solution_" + DSC::toString(current_refinement_) + "_factor_0_" + DSC::toString(ii), false);
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
      assert(current_solution_vector_);
      compute_reference_solution();
      assert(reference_discretization_);
      // compute error
      if (test_case_.provides_exact_solution()) {
        compute_exact_solution_vector();
        assert(exact_solution_vector_);
        std::unique_ptr< VectorType > difference = DSC::make_unique< VectorType >(*exact_solution_vector_);
        for (size_t ii = 0; ii < difference->size(); ++ii) {
          assert(DSC::FloatCmp::eq(difference->operator[](ii).first, current_solution_vector_->operator[](ii).first) && "Time steps must be the same");
          difference->operator[](ii).second = difference->operator[](ii).second - current_solution_vector_->operator[](ii).second;
        }
        return compute_norm(test_case_.reference_grid_view(), *difference, type);
      } else {
        assert(reference_solution_vector_);
        std::unique_ptr< VectorType > difference = DSC::make_unique< VectorType >(*reference_solution_vector_);
        for (size_t ii = 0; ii < difference->size(); ++ii) {
          assert(DSC::FloatCmp::eq(difference->operator[](ii).first, current_solution_vector_->operator[](ii).first) && "Time steps must be the same");
          difference->operator[](ii).second = difference->operator[](ii).second - current_solution_vector_->operator[](ii).second;
        }
        return compute_norm(test_case_.reference_grid_view(), *difference, type);
      }
    } else {
      assert(current_solution_vector_on_level_);
      return estimate(*current_solution_vector_on_level_, type);
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
                                                                                     test_case_.reference_level()));
      reference_solution_vector_ = Stuff::Common::make_unique< VectorType >(reference_discretization_->solve(test_case_.problem().is_linear()));
      reference_solution_computed_ = true;
      if (true) {
        DiscreteFunctionType tmp_discrete_func(reference_discretization_->fv_space(),
                                               reference_solution_vector_->operator[](0).second,
                                               "reference solution");
        for (size_t ii = 0; ii < reference_solution_vector_->size(); ++ii) {
          tmp_discrete_func.vector() = reference_solution_vector_->operator[](ii).second;
          tmp_discrete_func.template visualize_factor< 0 >(visualize_prefix_ + "_reference" + "_factor_0_" + DSC::toString(ii), false);
          tmp_discrete_func.template visualize_factor< 1 >(visualize_prefix_ + "_reference" + "_factor_1_" + DSC::toString(ii), false);
          tmp_discrete_func.template visualize_factor< 2 >(visualize_prefix_ + "_reference" + "_factor_2_" + DSC::toString(ii), false);
        }
      }
    }
  } // ... compute_reference_solution()

  void compute_exact_solution_vector()
  {
    if (!exact_solution_vector_computed_) {
      compute_reference_solution();
      const auto exact_solution = test_case_.exact_solution();
      VectorType discrete_exact_solution;
      for (size_t ii = 0; ii < reference_solution_vector_->size(); ++ii) {
        const double time = reference_solution_vector_->operator[](ii).first;
        const auto exact_solution_at_time = exact_solution->evaluate_at_time(time);
        DiscreteFunctionType discrete_exact_solution_at_time(reference_discretization_->fv_space(), "exact solution at time " + DSC::toString(time));
        project(*exact_solution_at_time, discrete_exact_solution_at_time);
        discrete_exact_solution.emplace_back(std::make_pair(time, discrete_exact_solution_at_time.vector()));
      }
      if (true) {
        DiscreteFunctionType tmp_discrete_func(reference_discretization_->fv_space(),
                                               discrete_exact_solution[0].second,
                                               "exact solution");
        for (size_t ii = 0; ii < discrete_exact_solution.size(); ++ii) {
          tmp_discrete_func.vector() = discrete_exact_solution[ii].second;
          tmp_discrete_func.template visualize_factor< 0 >(visualize_prefix_ + "_exact_solution" + "_factor_0_" + DSC::toString(ii), false);
          tmp_discrete_func.template visualize_factor< 1 >(visualize_prefix_ + "_exact_solution" + "_factor_1_" + DSC::toString(ii), false);
          tmp_discrete_func.template visualize_factor< 2 >(visualize_prefix_ + "_exact_solution" + "_factor_2_" + DSC::toString(ii), false);
        }
      }
      exact_solution_vector_ = Stuff::Common::make_unique< VectorType >(discrete_exact_solution);
//      std::cout << "norm: " << compute_norm(test_case_.reference_grid_view(), *exact_solution_vector_, "L1") << std::endl;
      exact_solution_vector_computed_ = true;
    }
  }

  bool is_norm(const std::string type) const
  {
    const auto norms = available_norms();
    return std::find(norms.begin(), norms.end(), type) != norms.end();
  }

  virtual std::vector< std::string > available_norms() const = 0;

  virtual std::vector< std::string > available_estimators() const = 0;

  virtual double estimate(const VectorType& vector, const std::string type) const = 0;

  virtual double compute_norm(const GridViewType& grid_view,
                              const VectorType& function,
                              const std::string type) const = 0;

  TestCaseType& test_case_;
  size_t current_refinement_;
  size_t last_computed_refinement_;
  mutable std::vector< double > grid_widths_;
  double time_to_solution_;
  bool reference_solution_computed_;
  bool exact_solution_vector_computed_;
  std::unique_ptr< DiscretizationType > current_discretization_;
  std::unique_ptr< VectorType > current_solution_vector_on_level_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
  std::unique_ptr< VectorType > reference_solution_vector_;
  std::unique_ptr< VectorType > current_solution_vector_;
  std::unique_ptr< VectorType > exact_solution_vector_;
  const std::string visualize_prefix_;
}; // class NonStationaryEocStudy


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_NONSTATIONARY_EOCSTUDY_HH
