// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH
#define DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH

#include <atomic>
#include <cmath>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/filters/intersection.hh>
#include <dune/xt/grid/integrals.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/functions/derivatives.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/laplace-ipdg.hh>
#include <dune/gdt/local/integrands/ipdg.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/norms.hh>
#include <dune/gdt/operators/constant.hh>
#include <dune/gdt/operators/laplace-ipdg-flux-reconstruction.hh>
#include <dune/gdt/operators/lincomb.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/oswald-interpolation.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \todo Add treatment of nonzero Dirichlet boundary values
 * \todo Add treatment of Neumann boundary values
 */
template <class G, XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class StationaryDiffusionIpdgEocStudy
  : public StationaryEocStudy<typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type,
                              1,
                              la>
{
protected:
  using BaseType =
      StationaryEocStudy<typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type, 1, la>;

  static const size_t d = BaseType::d;
  using typename BaseType::DF;
  using typename BaseType::GP;
  using typename BaseType::GV;
  using typename BaseType::I;
  using typename BaseType::M;
  using typename BaseType::O;
  using typename BaseType::S;
  using typename BaseType::V;

public:
  using typename BaseType::E;

  StationaryDiffusionIpdgEocStudy(const double& symmetry_prefactor,
                                  const double& inner_penalty,
                                  const double& dirichlet_penalty,
                                  const std::function<double(const I&)>& intersection_diameter =
                                      LocalIPDGIntegrands::internal::default_intersection_diameter<I>())
    : BaseType()
    , space_type_("")
    , symmetry_prefactor_(symmetry_prefactor)
    , inner_penalty_(inner_penalty)
    , dirichlet_penalty_(dirichlet_penalty)
    , intersection_diameter_(intersection_diameter)
  {}

protected:
  using FF = XT::Functions::GridFunctionInterface<E>;
  using FT = XT::Functions::GridFunctionInterface<E, d, d>;

  virtual const XT::Grid::BoundaryInfo<I>& boundary_info() const = 0;

  virtual const FT& diffusion() const = 0;

  virtual const FF& force() const = 0;

  virtual const FT& weight_function() const = 0;

  std::vector<std::string> norms() const override
  {
    auto nrms = BaseType::norms();
    nrms.push_back("eta_NC");
    nrms.push_back("eta_R");
    nrms.push_back("eta_DF");
    return nrms;
  }

  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level,
          const std::vector<std::string>& actual_norms,
          const std::vector<std::pair<std::string, std::string>>& actual_estimates,
          const std::vector<std::string>& actual_quantities) override
  {
    auto& self = *this;
    // compute the quantities/norms/estimates we known about and remove them from the todos
    // store the data in current_data_, BaseType::compute will return that
    auto remaining_norms = actual_norms;
    auto remaining_estimates = actual_estimates;
    auto remaining_quantities = actual_quantities;
    if (self.current_refinement_ != refinement_level)
      self.discretization_info(refinement_level);
    DUNE_THROW_IF(!self.current_space_, InvalidStateException, "");
    // compute current solution
    const auto& current_space = *self.current_space_;
    // visualize
    if (DXTC_TEST_CONFIG_GET("setup.visualize", false)) {
      const std::string prefix = XT::Common::Test::get_unique_test_name() + "_problem_";
      const std::string postfix = "_ref_" + XT::Common::to_string(refinement_level);
      //      self.diffusion_factor().visualize(current_space.grid_view(), prefix + "diffusion_factor" + postfix);
      self.diffusion().visualize(current_space.grid_view(), prefix + "diffusion" + postfix);
      self.force().visualize(current_space.grid_view(), prefix + "force" + postfix);
      //      self.dirichlet().visualize(current_space.grid_view(), prefix + "dirichlet" + postfix);
      //      self.neumann().visualize(current_space.grid_view(), prefix + "neumann" + postfix);
    }
    Timer timer;
    const auto solution = make_discrete_function(current_space, self.solve(current_space));
    // only set time if this did not happen in solve()
    if (self.current_data_["quantity"].count("time to solution (s)") == 0)
      self.current_data_["quantity"]["time to solution (s)"] = timer.elapsed();
    for (auto norm_it = remaining_norms.begin(); norm_it != remaining_norms.end(); /*Do not increment here ...*/) {
      const auto norm_id = *norm_it;
      if (norm_id == "eta_NC") {
        norm_it = remaining_norms.erase(norm_it); // ... but rather here ...
        // compute estimate
        auto oswald_interpolation_operator = make_oswald_interpolation_operator<M>(
            current_space.grid_view(), current_space, current_space, boundary_info());
        oswald_interpolation_operator.assemble(/*parallel=*/true);
        const auto h1_interpolation = oswald_interpolation_operator.apply(solution);
        self.current_data_["norm"][norm_id] =
            laplace_norm(current_space.grid_view(), /*weight=*/diffusion(), solution - h1_interpolation);
      } else if (norm_id == "eta_R") {
        norm_it = remaining_norms.erase(norm_it); // ... or here ...
        // compute estimate
        auto rt_space = make_raviart_thomas_space(current_space.grid_view(), current_space.max_polorder() - 1);
        auto reconstruction_op = make_laplace_ipdg_flux_reconstruction_operator<M>(current_space.grid_view(),
                                                                                   current_space,
                                                                                   rt_space,
                                                                                   symmetry_prefactor_,
                                                                                   inner_penalty_,
                                                                                   dirichlet_penalty_,
                                                                                   this->diffusion(),
                                                                                   this->weight_function(),
                                                                                   intersection_diameter_);
        auto flux_reconstruction = reconstruction_op.apply(solution);
        double eta_R_2 = 0.;
        std::mutex eta_R_2_mutex;
        auto walker = XT::Grid::make_walker(current_space.grid_view());
        walker.append(
            []() {},
            [&](const auto& element) {
              auto local_df = this->diffusion().local_function();
              local_df->bind(element);
              auto local_force = this->force().local_function();
              local_force->bind(element);
              auto local_flux = flux_reconstruction.local_function();
              local_flux->bind(element);
              auto flux_divergence = XT::Functions::divergence(*local_flux);
              flux_divergence.bind(element);
              // approximate minimum eigenvalue of the diffusion over the element ...
              double min_EV = std::numeric_limits<double>::max();
              // ... which we do by evaluating at some quadrature points
              for (auto&& quadrature_point : QuadratureRules<double, d>::rule(element.type(), local_df->order() + 3)) {
                auto diff = local_df->evaluate(quadrature_point.position());
                auto eigen_solver =
                    XT::LA::make_eigen_solver(diff,
                                              {{"type", XT::LA::EigenSolverOptions<decltype(diff)>::types().at(0)},
                                               {"assert_positive_eigenvalues", "1e-15"}});
                min_EV = std::min(min_EV, eigen_solver.min_eigenvalues(1).at(0));
              }
              DUNE_THROW_IF(!(min_EV > 0.),
                            Exceptions::integrand_error,
                            "The minimum eigenvalue of a positiv definite matrix must not be negative!"
                                << "\n\nmin_EV = " << min_EV);
              auto L2_norm_2 = LocalElementIntegralBilinearForm<E>(LocalProductIntegrand<E>(),
                                                                   /*over_integrate=*/3)
                                   .apply2(*local_force - flux_divergence, *local_force - flux_divergence)[0][0];
              const auto h = XT::Grid::diameter(element);
              const auto C_P = 1. / (M_PIl * M_PIl); // Poincare constant (known for simplices/cubes)
              std::lock_guard<std::mutex> lock(eta_R_2_mutex);
              eta_R_2 += (C_P * h * h * L2_norm_2) / min_EV;
            },
            []() {});
        walker.walk(/*parallel=*/true);
        self.current_data_["norm"][norm_id] = std::sqrt(eta_R_2);
      } else if (norm_id == "eta_DF") {
        norm_it = remaining_norms.erase(norm_it); // ... or here ...
        // compute estimate
        auto rt_space = make_raviart_thomas_space(current_space.grid_view(), current_space.max_polorder() - 1);
        auto reconstruction_op = make_laplace_ipdg_flux_reconstruction_operator<M>(current_space.grid_view(),
                                                                                   current_space,
                                                                                   rt_space,
                                                                                   symmetry_prefactor_,
                                                                                   inner_penalty_,
                                                                                   dirichlet_penalty_,
                                                                                   this->diffusion(),
                                                                                   this->weight_function(),
                                                                                   intersection_diameter_);
        auto flux_reconstruction = reconstruction_op.apply(solution);
        double eta_DF_2 = 0.;
        std::mutex eta_DF_2_mutex;
        auto walker = XT::Grid::make_walker(current_space.grid_view());
        walker.append(
            []() {},
            [&](const auto& element) {
              auto local_df = this->diffusion().local_function();
              local_df->bind(element);
              auto local_solution = solution.local_function();
              local_solution->bind(element);
              auto local_reconstruction = flux_reconstruction.local_function();
              local_reconstruction->bind(element);
              auto result = XT::Grid::element_integral(
                  element,
                  [&](const auto& xx) {
                    const auto diff = local_df->evaluate(xx);
                    const auto diff_inv = XT::LA::invert_matrix(diff);
                    const auto solution_grad = local_solution->jacobian(xx)[0];
                    const auto flux_rec = local_reconstruction->evaluate(xx);
                    auto difference = diff * solution_grad + flux_rec;
                    return (diff_inv * difference) * difference;
                  },
                  std::max(local_df->order() + std::max(local_solution->order() - 1, 0), local_reconstruction->order())
                      + /*over_integrate=*/3);
              std::lock_guard<std::mutex> lock(eta_DF_2_mutex);
              eta_DF_2 += result;
            },
            []() {});
        walker.walk(/*parallel=*/true);
        self.current_data_["norm"][norm_id] = std::sqrt(eta_DF_2);
      } else
        ++norm_it; // ... or finally here.
    } // norms
    // let the Base compute the rest
    return BaseType::compute(refinement_level, remaining_norms, remaining_estimates, remaining_quantities);
  } // ... compute(...)

  std::unique_ptr<S> make_space(const GP& current_grid) override
  {
    if (space_type_ == "fv")
      return std::make_unique<FiniteVolumeSpace<GV>>(current_grid.leaf_view());
    else if (space_type_.size() >= 4 && space_type_.substr(0, 4) == "dg_p") {
      const auto order = XT::Common::from_string<int>(space_type_.substr(4));
      return std::make_unique<DiscontinuousLagrangeSpace<GV>>(current_grid.leaf_view(), order);
    } else if (space_type_.size() >= 4 && space_type_.substr(0, 4) == "cg_p") {
      const auto order = XT::Common::from_string<int>(space_type_.substr(4));
      return std::make_unique<ContinuousLagrangeSpace<GV>>(current_grid.leaf_view(), order);
    } else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "space_type_ = " << space_type_);
      return nullptr;
    }
  } // ... make_space(...)

  std::unique_ptr<O> make_residual_operator(const S& space) override
  {
    // define lhs operator (has to be a pointer to allow the residual operator to manage the memory in the end)
    auto lhs_op = std::make_unique<MatrixOperator<M, GV>>(make_matrix_operator<M>(
        space,
        (space_type_.size() >= 2 && space_type_.substr(0, 2) == "cg") ? Stencil::element
                                                                      : Stencil::element_and_intersection));
    // - volume term
    lhs_op->append(LocalElementIntegralBilinearForm<E>(LocalLaplaceIntegrand<E>(this->diffusion())));
    // - inner faces
    lhs_op->append(LocalCouplingIntersectionIntegralBilinearForm<I>(
                       LocalLaplaceIPDGIntegrands::InnerCoupling<I>(1., this->diffusion(), this->weight_function())),
                   {},
                   XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
    lhs_op->append(LocalCouplingIntersectionIntegralBilinearForm<I>(LocalIPDGIntegrands::InnerPenalty<I>(
                       inner_penalty_, this->weight_function(), intersection_diameter_)),
                   {},
                   XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
    // - Dirichlet faces
    lhs_op->append(
        LocalIntersectionIntegralBilinearForm<I>(
            LocalLaplaceIPDGIntegrands::DirichletCoupling<I>(1., this->diffusion())),
        {},
        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(this->boundary_info(), new XT::Grid::DirichletBoundary()));
    lhs_op->append(
        LocalIntersectionIntegralBilinearForm<I>(LocalIPDGIntegrands::BoundaryPenalty<I>(
            dirichlet_penalty_, this->weight_function(), intersection_diameter_)),
        {},
        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(this->boundary_info(), new XT::Grid::DirichletBoundary()));
    // define rhs functional
    auto rhs_func = make_vector_functional<V>(space);
    rhs_func.append(LocalElementIntegralFunctional<E>(LocalProductIntegrand<E>().with_ansatz(this->force())));
    // ... add Dirichlet here
    // (if we add something here, the oswald interpolation in compute() needs to be adapted accordingly!)
    // ... add Neumann here
    // assemble everything in one grid walk
    lhs_op->append(rhs_func);
    lhs_op->assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    // build residual operator
    auto residual_op = std::make_unique<ConstLincombOperator<M, GV>>(space, space);
    residual_op->add(lhs_op.release(), 1.);
    residual_op->add(new ConstantOperator<M, GV>(space, space, new V(std::move(rhs_func.vector()))), -1);
    return residual_op;
  } // ... make_residual_operator(...)

  std::string space_type_;
  const double symmetry_prefactor_;
  const double inner_penalty_;
  const double dirichlet_penalty_;
  const std::function<double(const I&)> intersection_diameter_;
}; // class StationaryDiffusionIpdgEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_DIFFUSION_IPDG_HH
