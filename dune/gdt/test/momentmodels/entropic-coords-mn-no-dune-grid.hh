// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_NO_DUNE_GRID_HH
#define DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_NO_DUNE_GRID_HH

#include <chrono>
#include <numeric>

#include <boost/align/aligned_allocator.hpp>
#include <boost/align/is_aligned.hpp>
#include <boost/align/assume_aligned.hpp>

#include <dune/xt/common/numeric.hh>
#include <dune/xt/common/parallel/threadmanager.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/numerical-fluxes/kinetic.hh>
#include <dune/gdt/local/operators/advection-fv.hh>
#include <dune/gdt/local/operators/generic.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/operators/advection-fv-entropybased.hh>
#include <dune/gdt/operators/localizable-operator.hh>
#include <dune/gdt/operators/reconstruction/linear_kinetic.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/entropysolver.hh>
#include <dune/gdt/test/momentmodels/hessianinverter.hh>
#include <dune/gdt/test/momentmodels/density_evaluations.hh>
#include <dune/gdt/tools/timestepper/entropic-timestepper.hh>
#include <dune/gdt/tools/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/tools/timestepper/fractional-step.hh>
#include <dune/gdt/tools/timestepper/matrix-exponential-kinetic-isotropic.hh>

#include <dune/gdt/test/momentmodels/kineticequation.hh>

#include "pn-discretization.hh"

template <class DiscreteFunctionType>
struct ValuesVector
{
  static constexpr size_t dimRange = DiscreteFunctionType::r;
  using ThisType = ValuesVector;
  using RangeType = Dune::FieldVector<double, dimRange>;

  ValuesVector(const DiscreteFunctionType& discrete_function)
    : discrete_function_(discrete_function)
    , values_(discrete_function_.space().grid_view().size(0), RangeType(0))
  {
    const auto& fv_space = discrete_function_.space();
    const auto& mapper = fv_space.mapper();
    const auto& grid_view = fv_space.grid_view();
    // TODO: Only works if YaspGrid numbers the entities in the same way as we do. Check if this is true for two or more
    // dimensions.
    for (auto&& entity : Dune::elements(grid_view)) {
      for (size_t ii = 0; ii < dimRange; ++ii)
        values_[grid_view.indexSet().index(entity)][ii] =
            discrete_function_.dofs().vector()[mapper.global_index(entity, ii)];
    }
  }

  ThisType& operator+=(const ThisType& other)
  {
    assert(other.size() == size());
    for (size_t ii = 0; ii < size(); ++ii)
      values_[ii] += other.values_[ii];
    return *this;
  }

  ThisType& operator-=(const ThisType& other)
  {
    assert(other.size() == size());
    for (size_t ii = 0; ii < size(); ++ii)
      values_[ii] -= other.values_[ii];
    return *this;
  }

  ThisType& operator*=(const double val)
  {
    double* entries = data();
    const size_t num_entries = num_elements();
    for (size_t ii = 0; ii < num_entries; ++ii)
      entries[ii] *= val;
    return *this;
  }

  void set_zero()
  {
    double* entries = data();
    BOOST_ALIGN_ASSUME_ALIGNED(entries, 64);
    std::fill(entries, entries + num_elements(), 0.);
  }

  void copy_values(const ThisType& other)
  {
    double* entries = data();
    const double* other_entries = other.data();
    BOOST_ALIGN_ASSUME_ALIGNED(entries, 64);
    BOOST_ALIGN_ASSUME_ALIGNED(other_entries, 64);
    const size_t num_entries = num_elements();
    for (size_t ii = 0; ii < num_entries; ++ii)
      entries[ii] = other_entries[ii];
    // std::copy(other_entries, other_entries + num_entries, entries);
  }

  void axpy(const double val, const ThisType& other)
  {
    assert(other.size() == size());
    double* entries = data();
    const double* other_entries = other.data();
    BOOST_ALIGN_ASSUME_ALIGNED(entries, 64);
    BOOST_ALIGN_ASSUME_ALIGNED(other_entries, 64);
    const size_t num_entries = num_elements();
    for (size_t ii = 0; ii < num_entries; ++ii)
      entries[ii] += other_entries[ii] * val;
  }

  RangeType& operator[](const size_t ii)
  {
    return values_[ii];
  }

  const RangeType& operator[](const size_t ii) const
  {
    return values_[ii];
  }

  size_t size() const
  {
    return values_.size();
  }

  size_t num_elements() const
  {
    return size() * dimRange;
  }

  double* data()
  {
    return &(values_[0][0]);
  }

  const double* data() const
  {
    return &(values_[0][0]);
  }

  std::vector<RangeType, boost::alignment::aligned_allocator<RangeType, 64>>& values()
  {
    return values_;
  }

  const std::vector<RangeType, boost::alignment::aligned_allocator<RangeType, 64>>& values() const
  {
    return values_;
  }

  // TODO: Only works if YaspGrid numbers the entities in the same way as we do. Check if this is true for two or more
  // dimensions.
  template <class VisualizerType>
  void visualize(const std::string& filename, const VisualizerType& visualizer) const
  {
    const auto& fv_space = discrete_function_.space();
    const auto& mapper = fv_space.mapper();
    const auto& grid_view = fv_space.grid_view();
    for (auto&& entity : Dune::elements(grid_view)) {
      for (size_t ii = 0; ii < dimRange; ++ii)
        discrete_function_.dofs().vector()[mapper.global_index(entity, ii)] =
            values_[grid_view.indexSet().index(entity)][ii];
    }
    discrete_function_.visualize(grid_view, filename, false, Dune::VTK::appendedraw, {}, visualizer);
  }

  mutable DiscreteFunctionType discrete_function_;
  std::vector<RangeType, boost::alignment::aligned_allocator<RangeType, 64>> values_;
};

template <size_t d, size_t r, class ProblemType, class EntropyFluxType, class DiscreteFunctionType>
struct EntropicMnOperator
{
  using DomainType = Dune::FieldVector<double, d>;
  using RangeType = Dune::FieldVector<double, r>;
  using ValuesType = ValuesVector<DiscreteFunctionType>;
  using GridViewType = typename DiscreteFunctionType::SpaceType::GridViewType;

  EntropicMnOperator(const GridViewType& grid_view,
                     const Dune::FieldVector<size_t, d> grid_sizes,
                     const DomainType lower_left,
                     const DomainType upper_right,
                     const ProblemType& problem,
                     EntropyFluxType& entropy_flux,
                     const double psi_min)
    : grid_size_(
          std::accumulate(grid_sizes.begin(), grid_sizes.end(), 1., [](const auto& a, const auto& b) { return a * b; }))
    , grid_(grid_size_)
    , sigma_a_(grid_size_)
    , sigma_s_(grid_size_)
    , Q_(grid_size_)
    , lower_left_(lower_left)
    , upper_right_(upper_right)
    , problem_(problem)
    , entropy_flux_(entropy_flux)
    , psi_min_(psi_min)
    , reg_indicators_(grid_size_)
    , h_((upper_right_[0] - lower_left_[0]) / grid_sizes[0])
    , u_iso_(entropy_flux_.basis_functions().u_iso())
    , basis_integrated_(entropy_flux_.basis_functions().integrated())
    , alpha_one_(entropy_flux_.basis_functions().alpha_one())
    , left_boundary_flux_(problem.kinetic_boundary_flux(lower_left_, -1., 0))
    , right_boundary_flux_(problem.kinetic_boundary_flux(upper_right_, 1., 0))
  {
    // assume equally-spaced grid for now
    for (size_t dd = 1; dd < d; ++dd)
      assert(Dune::XT::Common::FloatCmp::eq((upper_right_[dd] - lower_left_[dd]) / grid_sizes[dd], h_));
    // store grid centers
    if (d == 1) {
      for (size_t nn = 0; nn < grid_size_; ++nn)
        grid_[nn] = lower_left_[0] + (nn + 0.5) * h_;
    } else {
      DUNE_THROW(Dune::NotImplemented, "");
    }
    // evaluate and store problem parameters
    const auto sigma_a_func = problem_.sigma_a();
    const auto sigma_s_func = problem_.sigma_s();
    const auto Q_func = problem_.Q();
    for (size_t nn = 0; nn < grid_size_; ++nn) {
      sigma_a_[nn] = sigma_a_func->evaluate(grid_[nn])[0];
      sigma_s_[nn] = sigma_s_func->evaluate(grid_[nn])[0];
      Q_[nn] = Q_func->evaluate(grid_[nn])[0];
    }
    entropy_flux_.prepare_storage(grid_view);
    // precompute evaluations of boundary distribution at quadrature points
    entropy_flux_.store_boundary_evaluations(
        problem_.boundary_distribution()(grid_[0] - h_ / 2), /*entity_index*/ 0, /*intersection.indexInInside()*/ 0);
    entropy_flux_.store_boundary_evaluations(problem_.boundary_distribution()(grid_.back() + DomainType(h_ / 2)),
                                             /*entity_index*/ grid_.size() - 1,
                                             /*intersection.indexInInside()*/ 1);
  }

  void apply(ValuesType& source_alpha, ValuesType& range_alpha, const bool regularize)
  {
    apply_density_op(source_alpha);
    ValuesType& update_u = range_alpha;
    update_u.set_zero();
    apply_advection_op(update_u);
    apply_rhs_op(update_u);
    std::fill(reg_indicators_.begin(), reg_indicators_.end(), false);
    // range_alpha contains the u updates now, transform to alpha updates next by applying the inverse hessian
    apply_inverse_hessian(range_alpha, regularize);
  }

  void apply_density_op(ValuesType& source_alpha)
  {
    for (size_t nn = 0; nn < grid_.size(); ++nn) {
      entropy_flux_.store_evaluations(nn, source_alpha[nn], psi_min_, true);
    }
    entropy_flux_.set_eta_ast_pointers();
  }

  // TODO: only works for one dimension by now
  // No need for source values as we use the precomputed density values from entropy_flux_
  void apply_advection_op(ValuesType& range_u)
  {
    const double h_inv = 1. / h_;
    const auto& densities = entropy_flux_.exp_evaluations();
    const auto& boundary_densities = entropy_flux_.boundary_distribution_evaluations();
    // left boundary
    range_u[0].axpy(h_inv, left_boundary_flux_);
    // loop over entities
    const size_t grid_size = grid_.size();
    for (size_t nn = 0; nn < grid_size; ++nn) {
      entropy_flux_.apply_kinetic_flux_with_kinetic_reconstruction(
          h_inv,
          nn == 0 ? boundary_densities[0][0][0] : densities[nn - 1],
          densities[nn],
          nn == grid_size - 1 ? boundary_densities[grid_size - 1][0][1] : densities[nn + 1],
          nn == 0 ? nullptr : &(range_u[nn - 1]),
          &(range_u[nn]),
          nn == grid_size - 1 ? nullptr : &(range_u[nn + 1]),
          0);
    }
    // right boundary
    range_u[grid_size - 1].axpy(h_inv, right_boundary_flux_);
  }

  // No need for source values as we use the precomputed density values from entropy_flux_
  void apply_rhs_op(ValuesType& range_u)
  {
    RangeType u_elem;
    for (size_t nn = 0; nn < grid_.size(); ++nn) {
      entropy_flux_.get_u(nn, u_elem);
      range_u[nn].axpy(-(sigma_a_[nn] + sigma_s_[nn]), u_elem);
      range_u[nn].axpy(entropy_flux_.basis_functions().density(u_elem) * sigma_s_[nn], u_iso_);
      range_u[nn].axpy(Q_[nn], basis_integrated_);
    };
  }

  void apply_inverse_hessian(ValuesType& update_u, const bool regularize)
  {
    for (size_t nn = 0; nn < grid_size_; ++nn) {
      try {
        entropy_flux_.apply_inverse_hessian(nn, update_u[nn]);
        for (auto&& entry : update_u[nn])
          if (std::isnan(entry) || std::isinf(entry)) {
            DUNE_THROW(Dune::MathError, "Hessian");
          }
      } catch (const Dune::MathError& e) {
        if (regularize) {
          reg_indicators_[nn] = true;
          return;
        } else
          throw e;
      }
    } // nn
  }

  const std::vector<bool>& reg_indicators() const
  {
    return reg_indicators_;
  }

  void regularize(const double rr, ValuesType& alpha)
  {
    const auto& basis_functions = entropy_flux_.basis_functions();
    for (size_t nn = 0; nn < grid_size_; ++nn) {
      if (reg_indicators_[nn]) {
        std::cout << "Regularized on entity " << nn << " with r = " << rr << std::endl;
        const double old_density = basis_functions.density(entropy_flux_.get_u(alpha[nn]));
        const RangeType alpha_iso = basis_functions.alpha_iso(old_density);
        alpha[nn] *= 1 - rr;
        alpha[nn].axpy(rr, alpha_iso);
        const double reg_density = basis_functions.density(entropy_flux_.get_u(alpha[nn]));
        const double factor = std::log(old_density / reg_density);
        alpha[nn].axpy(factor, alpha_one_);
      }
    }
  }

  double psi_min() const
  {
    return psi_min_;
  }

  const EntropyFluxType& entropy_flux() const
  {
    return entropy_flux_;
  }

  const size_t grid_size_;
  std::vector<DomainType> grid_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_s_;
  std::vector<double> Q_;
  DomainType lower_left_;
  DomainType upper_right_;
  const ProblemType& problem_;
  EntropyFluxType& entropy_flux_;
  const double psi_min_;
  std::vector<bool> reg_indicators_;
  const double h_;
  const RangeType u_iso_;
  const RangeType basis_integrated_;
  const RangeType alpha_one_;
  const RangeType left_boundary_flux_;
  const RangeType right_boundary_flux_;
};

template <class TestCaseType>
struct HyperbolicEntropicCoordsMnNoDuneGridDiscretization
{
  // returns: (l1norm, l2norm, linfnorm, MPI rank)
  static std::pair<Dune::FieldVector<double, 3>, int> run(size_t num_save_steps = 1,
                                                          size_t num_output_steps = 0,
                                                          size_t quad_order = size_t(-1),
                                                          size_t quad_refinements = size_t(-1),
                                                          std::string grid_size = "",
                                                          size_t overlap_size = 2,
                                                          double t_end = 0.,
                                                          std::string filename = "",
                                                          bool /*disable_thread_cache*/ = false)
  {
    using namespace Dune;
    using namespace Dune::GDT;

    //******************* get typedefs and constants from ProblemType **********************//
    using MomentBasis = typename TestCaseType::MomentBasis;
    using DiscreteFunctionType = typename TestCaseType::DiscreteFunctionType;
    using GridType = typename TestCaseType::GridType;
    using SpaceType = typename TestCaseType::SpaceType;
    using AdvectionSourceSpaceType = typename TestCaseType::AdvectionSourceSpaceType;
    using GV = typename TestCaseType::GridViewType;
    //    using E = XT::Grid::extract_entity_t<GV>;
    using ProblemType = typename TestCaseType::ProblemType;
    using RangeFieldType = typename MomentBasis::RangeFieldType;
    static constexpr size_t dimDomain = MomentBasis::dimDomain;
    static constexpr size_t dimRange = MomentBasis::dimRange;
    using DomainType = FieldVector<RangeFieldType, dimDomain>;

    //******************* create grid and FV space ***************************************
    auto grid_config = ProblemType::default_grid_cfg();
    if (!grid_size.empty())
      grid_config["num_elements"] = grid_size;
    grid_config["overlap_size"] = XT::Common::to_string(overlap_size);
    const auto grid_ptr =
        Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(grid_config, MPIHelper::getCommunicator()).grid_ptr();
    assert(grid_ptr->comm().size() == 1 || grid_ptr->overlapSize(0) > 0);
    const GV grid_view(grid_ptr->leafGridView());
    const SpaceType fv_space(grid_view);
    const AdvectionSourceSpaceType advection_source_space(grid_view);

    //******************* create EquationType object ***************************************
    std::shared_ptr<const MomentBasis> basis_functions = std::make_shared<const MomentBasis>(
        quad_order == size_t(-1) ? MomentBasis::default_quad_order() : quad_order,
        quad_refinements == size_t(-1) ? MomentBasis::default_quad_refinements() : quad_refinements);
    const std::unique_ptr<ProblemType> problem_ptr =
        XT::Common::make_unique<ProblemType>(*basis_functions, grid_view, grid_config);
    const auto& problem = *problem_ptr;
    const auto initial_values_u = problem.initial_values();
    const auto boundary_values_u = problem.boundary_values();
    const auto boundary_distribution = problem.boundary_distribution();

    constexpr SlopeLimiterType slope =
        TestCaseType::reconstruction ? SlopeLimiterType::minmod : SlopeLimiterType::no_slope;
    // TestCaseType::reconstruction ? SlopeLimiterType::superbee : SlopeLimiterType::no_slope;
    using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GV, MomentBasis, slope>;
    using OldEntropyFluxType = EntropyBasedFluxFunction<GV, MomentBasis>;
    auto flux = problem.flux();
    auto* entropy_flux = dynamic_cast<OldEntropyFluxType*>(flux.get());
    auto analytical_flux = std::make_unique<EntropyFluxType>(*entropy_flux);

    // ***************** project initial values to discrete function *********************
    // create a discrete function for the solution
    DiscreteFunctionType u(fv_space, "u_initial");
    DiscreteFunctionType alpha(fv_space, "alpha_initial");
    // project initial values
    default_interpolation(*initial_values_u, u, grid_view);

    const auto u_local_func = u.local_discrete_function();
    const auto alpha_local_func = alpha.local_discrete_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> u_local;
    for (auto&& element : Dune::elements(grid_view)) {
      u_local_func->bind(element);
      alpha_local_func->bind(element);
      for (size_t ii = 0; ii < dimRange; ++ii)
        u_local[ii] = u_local_func->dofs().get_entry(ii);
      const auto alpha_local = analytical_flux->get_alpha(u_local);
      for (size_t ii = 0; ii < dimRange; ++ii)
        alpha_local_func->dofs().set_entry(ii, alpha_local[ii]);
    }

    // ******************** choose flux and rhs operator and timestepper ******************************************

    // *************** Calculate dx and initial dt **************************************
    Dune::XT::Grid::Dimensions<GV> dimensions(grid_view);
    RangeFieldType dx = dimensions.entity_width.max();
    if (dimDomain == 2)
      dx /= std::sqrt(2);
    if (dimDomain == 3)
      dx /= std::sqrt(3);

    // *********************** create operators and timesteppers ************************************
    const double min_acceptable_density = problem.psi_vac() / 10;
    using DomainSizesType = FieldVector<size_t, dimDomain>;
    using MnOperatorType = EntropicMnOperator<dimDomain, dimRange, ProblemType, EntropyFluxType, DiscreteFunctionType>;
    MnOperatorType entropic_mn_operator(grid_view,
                                        XT::Common::from_string<DomainSizesType>(grid_config["num_elements"]),
                                        XT::Common::from_string<DomainType>(grid_config["lower_left"]),
                                        XT::Common::from_string<DomainType>(grid_config["upper_right"]),
                                        problem,
                                        *analytical_flux,
                                        min_acceptable_density);
    using ValuesType = ValuesVector<DiscreteFunctionType>;
    ValuesType alpha_vals(alpha);
    using TimeStepperType =
        EntropicTimeStepper<MnOperatorType, EntropyFluxType, ValuesType, TimeStepperMethods::bogacki_shampine>;
    // EntropicTimeStepper<MnOperatorType, EntropyFluxType, ValuesType, TimeStepperMethods::dormand_prince>;

    if (XT::Common::is_zero(t_end))
      t_end = problem.t_end();

    if (!filename.empty())
      filename += "_";
    if (TestCaseType::reconstruction && slope == SlopeLimiterType::minmod)
      filename += "minmod_";
    else if (TestCaseType::reconstruction && slope == SlopeLimiterType::superbee)
      filename += "superbee_";
    filename += ProblemType::static_id();
    filename += "_grid_" + grid_config["num_elements"];
    filename += "_tend_" + XT::Common::to_string(t_end);
    filename += "_quad_" + XT::Common::to_string(quad_order);
    filename += TestCaseType::reconstruction ? "_ord2" : "_ord1";
    filename += "_" + basis_functions->mn_name();

    // ******************************** do the time steps ***********************************************************
    TimeStepperType timestepper(entropic_mn_operator, *analytical_flux, alpha_vals, 1.);

    auto begin_time = std::chrono::steady_clock::now();
    auto visualizer = std::make_unique<XT::Functions::GenericVisualizer<dimRange, 1, double>>(
        1, [&basis_functions, &analytical_flux](const int /*comp*/, const auto& val) {
          return basis_functions->density(analytical_flux->get_u(val));
        });

    // The hessian has entries in the order of psi_min, the inverse thus scales with 1/psi_min, and thus the timestep
    // should be psi_min to get an update of order 1
    double initial_dt = min_acceptable_density;
    timestepper.solve(t_end, initial_dt, num_save_steps, num_output_steps, true, filename, *visualizer);
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_diff = end_time - begin_time;
    if (grid_view.comm().rank() == 0)
      std::cout << "Solving took: " << XT::Common::to_string(time_diff.count(), 15) << " s" << std::endl;

    // auto ret = std::make_pair(FieldVector<double, 3>(0.), int(0));
    // double& l1norm = ret.first[0];
    // double& l2norm = ret.first[1];
    // double& linfnorm = ret.first[2];
    // ret.second = grid_view.comm().rank();
    // const auto& current_sol = timestepper.current_solution();
    // const auto local_sol = current_sol.local_function();
    // for (const auto& entity : elements(grid_view, Dune::Partitions::interior)) {
    //   local_sol->bind(entity);
    //   const auto val = local_sol->evaluate(entity.geometry().local(entity.geometry().center()));
    //   RangeFieldType psi = basis_functions->density(val);
    //   l1norm += std::abs(psi) * entity.geometry().volume();
    //   l2norm += std::pow(psi, 2) * entity.geometry().volume();
    //   linfnorm = std::max(std::abs(psi), linfnorm);
    // }
    // l1norm = grid_view.comm().sum(l1norm);
    // l2norm = grid_view.comm().sum(l2norm);
    // linfnorm = grid_view.comm().max(linfnorm);
    // l2norm = std::sqrt(l2norm);
    auto ret = std::make_pair(FieldVector<double, 3>(0.), int(0));
    return ret;
  }
};

template <class TestCaseType>
struct HyperbolicEntropicCoordsMnNoDuneGridTest
  : public HyperbolicEntropicCoordsMnNoDuneGridDiscretization<TestCaseType>
  , public ::testing::Test
{
  void run()
  {
    auto norms = HyperbolicEntropicCoordsMnNoDuneGridDiscretization<TestCaseType>::run(
                     DXTC_CONFIG.get("num_save_steps", 10),
                     -1,
                     TestCaseType::quad_order,
                     TestCaseType::quad_refinements,
                     DXTC_CONFIG.get("grid_size", ""),
                     2,
                     DXTC_CONFIG.get("t_end", TestCaseType::t_end),
                     "test_nodune",
                     Dune::GDT::is_full_moment_basis<typename TestCaseType::MomentBasis>::value)
                     .first;
    const double l1norm = norms[0];
    const double l2norm = norms[1];
    const double linfnorm = norms[2];
    using ResultsType = typename TestCaseType::ExpectedResultsType;
    EXPECT_NEAR(ResultsType::l1norm, l1norm, ResultsType::l1norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::l2norm, l2norm, ResultsType::l2norm * ResultsType::tol);
    EXPECT_NEAR(ResultsType::linfnorm, linfnorm, ResultsType::linfnorm * ResultsType::tol);
  }
};

#endif // DUNE_GDT_TEST_HYPERBOLIC_ENTROPIC_COORDS_MN_NO_DUNE_GRID_HH
