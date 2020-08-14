// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include "config.h"

#include <mutex>
#include <tuple>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/integrals.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/derivatives.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/spe10/model1.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/ipdg-flux-reconstruction.hh>
#include <dune/gdt/operators/oswald-interpolation.hh>
#include <dune/gdt/prolongations.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

using namespace Dune;
using namespace Dune::GDT;

#if HAVE_DUNE_ALUGRID
using G = ALU_2D_SIMPLEX_CONFORMING;
#elif HAVE_DUNE_UGGRID || HAVE_UG
using G = UG_2D;
#else
#  warning Falling back to cubic grid, results will not be reproduced but similar!
using G = YASP_2D_EQUIDISTANT_OFFSET;
#endif

using GP = XT::Grid::GridProvider<G>;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
static const constexpr size_t d = G::dimension;

using M = XT::LA::IstlRowMajorSparseMatrix<double>;
using V = XT::LA::IstlDenseVector<double>;

using DG = DiscontinuousLagrangeSpace<GV>;
using RTN = RaviartThomasSpace<GV>;
using ScalarDF = DiscreteFunction<V, GV>;
using VectorDF = DiscreteFunction<V, GV, d>;


std::unique_ptr<M> assemble_SWIPDG_matrix(const DG& space,
                                          const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                          const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                          const bool parallel)
{
  auto op = make_matrix_operator<M>(space, Stencil::element_and_intersection);
  op.append(LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(diffusion_factor, diffusion_tensor)));
  op.append(
      LocalCouplingIntersectionIntegralBilinearForm<I>(
          LocalEllipticIpdgIntegrands::Inner<I, double, LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
              diffusion_factor, diffusion_tensor)),
      {},
      XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
  op.append(LocalCouplingIntersectionIntegralBilinearForm<I>(
                LocalEllipticIpdgIntegrands::
                    DirichletBoundaryLhs<I, double, LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
                        diffusion_factor, diffusion_tensor)),
            {},
            XT::Grid::ApplyOn::BoundaryIntersections<GV>());
  op.assemble(parallel);
  return std::make_unique<M>(std::move(op.matrix()));
} // ... assemble_SWIPDG_matrix(...)


std::unique_ptr<V>
assemble_L2_vector(const DG& space, const XT::Functions::GridFunctionInterface<E>& force, const bool parallel)
{
  auto func = make_vector_functional<V>(space);
  func.append(LocalElementIntegralFunctional<E>(LocalElementProductIntegrand<E>().with_ansatz(force)));
  func.assemble(parallel);
  return std::make_unique<V>(std::move(func.vector()));
} // ... assemble_L2_vector(...)


std::unique_ptr<M>
assemble_energy_semi_product_matrix(const DG& space,
                                    const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                    const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                    const bool parallel)
{
  auto op = make_matrix_operator<M>(space, Stencil::element_and_intersection);
  op.append(LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(diffusion_factor, diffusion_tensor)));
  op.assemble(parallel);
  return std::make_unique<M>(std::move(op.matrix()));
} // ... assemble_energy_semi_product_matrix(...)


std::unique_ptr<M> assemble_DG_product_matrix(const DG& space,
                                              const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                              const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                              const bool parallel)
{
  auto op = make_matrix_operator<M>(space, Stencil::element_and_intersection);
  op.append(LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(diffusion_factor, diffusion_tensor)));
  op.append(LocalCouplingIntersectionIntegralBilinearForm<I>(
                LocalEllipticIpdgIntegrands::
                    InnerOnlyPenalty<I, double, LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
                        diffusion_factor, diffusion_tensor)),
            {},
            XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
  op.append(
      LocalCouplingIntersectionIntegralBilinearForm<I>(
          LocalEllipticIpdgIntegrands::
              DirichletBoundaryLhsOnlyPenalty<I, double, LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
                  diffusion_factor, diffusion_tensor)),
      {},
      XT::Grid::ApplyOn::BoundaryIntersections<GV>());
  op.assemble(parallel);
  return std::make_unique<M>(std::move(op.matrix()));
} // ... assemble_DG_product_matrix(...)


std::unique_ptr<V> compute_flux_reconstruction(const GP& grid,
                                               const DG& dg_space,
                                               const RTN& rtn_space,
                                               const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                               const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                               const V& dg_vec)
{
  auto op = make_ipdg_flux_reconstruction_operator<M, LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
      grid.leaf_view(), dg_space, rtn_space, diffusion_factor, diffusion_tensor);
  auto rtn_vec = op.apply(dg_vec);
  return std::make_unique<V>(std::move(rtn_vec));
}


double compute_local_conservation_error(const GP& grid,
                                        const VectorDF& flux,
                                        const XT::Functions::GridFunctionInterface<E>& rhs,
                                        const bool parallel)
{
  double error = 0.;
  std::mutex eror_mutex;
  auto walker = XT::Grid::make_walker(grid.leaf_view());
  walker.append([](/*prepare nothing*/) {},
                [&](const auto& element) {
                  auto local_flux = flux.local_function();
                  local_flux->bind(element);
                  auto local_rhs = rhs.local_function();
                  local_rhs->bind(element);
                  auto local_error = LocalElementIntegralFunctional<E>(
                                         [&](const auto&, const auto&) {
                                           return std::max(std::max(local_flux->order() - 1, 0), local_rhs->order());
                                         },
                                         [&](const auto&, const auto& xx, auto& result, const auto&) {
                                           auto flux_grads = local_flux->jacobian(xx);
                                           auto rhs_value = local_rhs->evaluate(xx);
                                           auto divergence = [](const auto& grad) {
                                             double div = 0.;
                                             for (size_t dd = 0; dd < d; ++dd)
                                               div += grad[dd][dd];
                                             return div;
                                           };
                                           result[0] = divergence(flux_grads) - rhs_value;
                                         },
                                         {})
                                         .apply(*local_rhs)[0];
                  std::lock_guard<std::mutex> lock(eror_mutex);
                  error += std::abs(local_error);
                },
                [](/*finalize nothing*/) {});
  walker.walk(parallel);
  return error;
} // ... compute_local_conservation_error(...)


std::tuple<double, double, double, double>
compute_estimate(const GP& grid,
                 const ScalarDF& pressure,
                 const VectorDF& flux,
                 const XT::Functions::GridFunctionInterface<E>& rhs,
                 const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                 const XT::Functions::GridFunctionInterface<E>& diffusion_factor_bar,
                 const XT::Functions::GridFunctionInterface<E>& diffusion_factor_hat,
                 const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                 const double& alpha_bar,
                 const double& alpha_hat,
                 const double& gamma_bar,
                 const int over_integrate,
                 const bool& parallel)
{
  auto diffusion = diffusion_factor * diffusion_tensor;
  auto diffusion_hat = diffusion_factor_hat * diffusion_tensor;
  double eta_2 = 0.;
  double eta_NC_2 = 0.;
  double eta_R_2 = 0.;
  double eta_DF_2 = 0.;
  std::mutex mutex;
  auto grid_view = grid.leaf_view();
  const XT::Grid::AllDirichletBoundaryInfo<I> boundary_info;
  const auto& dg_space = pressure.space();
  auto oswald_interpolation_operator =
      make_oswald_interpolation_operator<M>(grid_view, dg_space, dg_space, boundary_info);
  oswald_interpolation_operator.assemble(parallel);
  const auto conforming_pressure = oswald_interpolation_operator.apply(pressure);
  auto walker = XT::Grid::make_walker(grid.leaf_view());
  walker.append(
      [](/*prepare nothing*/) {},
      [&](const auto& element) {
        // prepare data functions
        auto p = pressure.local_function();
        auto cp = conforming_pressure.local_function();
        auto t = flux.local_function();
        auto div_t = XT::Functions::divergence(*t);
        auto f = rhs.local_function();
        auto df = diffusion.local_function();
        auto dfh = diffusion_hat.local_function();
        p->bind(element);
        cp->bind(element);
        t->bind(element);
        div_t.bind(element);
        f->bind(element);
        df->bind(element);
        dfh->bind(element);
        // eta_NC
        const double eta_NC_element_2 =
            LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(diffusion_factor_bar, diffusion_tensor),
                                                over_integrate)
                .apply2(*p - *cp, *p - *cp)[0][0];
        // eta_R
        // - approximate minimum eigenvalue of the diffusion over the element (evaluate at some points)
        double min_EV = std::numeric_limits<double>::max();
        for (auto&& quadrature_point :
             QuadratureRules<double, d>::rule(element.type(), dfh->order() + over_integrate)) {
          auto diff = dfh->evaluate(quadrature_point.position());
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
        const auto C_P = 1. / (M_PIl * M_PIl); // Poincare constant (known for simplices/cubes)
        const auto h = XT::Grid::diameter(element);
        auto L2_norm_2 = LocalElementIntegralBilinearForm<E>(LocalElementProductIntegrand<E>(), over_integrate)
                             .apply2(*f - div_t, *f - div_t)[0][0];
        const double eta_R_element_2 = (C_P * h * h * L2_norm_2) / min_EV;
        // eta_DF
        const double eta_DF_element_2 = XT::Grid::element_integral(
            element,
            [&](const auto& xx) {
              const auto diff = df->evaluate(xx);
              const auto diff_inv = XT::LA::invert_matrix(dfh->evaluate(xx));
              const auto pressure_grad = p->jacobian(xx)[0];
              const auto t_val = t->evaluate(xx);
              auto difference = diff * pressure_grad + t_val;
              return (diff_inv * difference) * difference;
            },
            std::max(df->order() + std::max(p->order() - 1, 0), t->order()) + over_integrate);
        // compute indicators and estimator
        std::lock_guard<std::mutex> lock(mutex);
        eta_NC_2 += eta_NC_element_2;
        eta_R_2 += eta_R_element_2;
        eta_DF_2 += eta_DF_element_2;
        eta_2 += gamma_bar * eta_NC_element_2
                 + std::pow(std::sqrt(eta_R_element_2) + std::sqrt(eta_DF_element_2) / std::sqrt(alpha_hat), 2);
      },
      [](/*finalize nothing*/) {});
  walker.walk(parallel);
  eta_2 /= alpha_bar;
  return std::make_tuple(std::sqrt(eta_2), std::sqrt(eta_NC_2), std::sqrt(eta_R_2), std::sqrt(eta_DF_2));
} // ... compute_estimate(...)


PYBIND11_MODULE(gamm_2019_talk_on_conservative_rb, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt.discretefunction");

  m.def(
      "visualize",
      [](GP& grid, XT::Functions::FunctionInterface<d>& func, const std::string& filename) {
        func.visualize(grid.leaf_view(), filename);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "scalar_function"_a,
      "filename"_a);
  m.def(
      "visualize",
      [](GP& grid, XT::Functions::FunctionInterface<d, d, d>& func, const std::string& filename) {
        func.visualize(grid.leaf_view(), filename);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "matrix_function"_a,
      "filename"_a);
  m.def(
      "visualize",
      [](GP& grid, XT::Functions::GridFunctionInterface<E>& func, const std::string& filename) {
        func.visualize(grid.leaf_view(), filename);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "scalar_function"_a,
      "filename"_a);
  m.def(
      "visualize",
      [](GP& grid, XT::Functions::GridFunctionInterface<E, d, d>& func, const std::string& filename) {
        func.visualize(grid.leaf_view(), filename);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "scalar_function"_a,
      "filename"_a);

  py::class_<GP> grid_provider(m, "GridProvider", "GridProvider");
  grid_provider.def(py::init([](const FieldVector<double, d>& lower_left,
                                const FieldVector<double, d>& upper_right,
                                const std::array<unsigned int, d>& num_elements) {
                      return new GP(XT::Grid::make_cube_grid<G>(lower_left, upper_right, num_elements));
                    }),
                    "lower_left"_a,
                    "upper_right"_a,
                    "num_elements"_a);
  grid_provider.def_property_readonly("num_elements", [](GP& self) { return self.leaf_view().indexSet().size(0); });
  grid_provider.def("refine", [](GP& self, const int num_refinements) { self.global_refine(num_refinements); });

  py::class_<DG> dg_space(m, "DiscontinuousLagrangeSpace", "DiscontinuousLagrangeSpace");
  dg_space.def(py::init([](GP& grid_provider, const int order) { return new DG(grid_provider.leaf_view(), order); }),
               "grid_provider"_a,
               "order"_a = 1);
  dg_space.def_property_readonly("dimDomain", [](DG& /*self*/) { return d; });
  dg_space.def_property_readonly("num_DoFs", [](DG& self) { return self.mapper().size(); });

  py::class_<RTN> rtn_space(m, "RaviartThomasSpace", "RaviartThomasSpace");
  rtn_space.def(py::init([](GP& grid_provider, const int order) { return new RTN(grid_provider.leaf_view(), order); }),
                "grid_provider"_a,
                "order"_a = 1);
  rtn_space.def_property_readonly("dimDomain", [](RTN& /*self*/) { return d; });
  rtn_space.def_property_readonly("num_DoFs", [](RTN& self) { return self.mapper().size(); });

  m.def(
      "make_discrete_function",
      [](DG& dg_space, V& vec, const std::string& name) { return ScalarDF(dg_space, vec, name); },
      "dg_space"_a,
      "DoF_vector"_a,
      "name"_a);
  m.def(
      "make_discrete_function",
      [](RTN& rtn_space, V& vec, const std::string& name) { return VectorDF(rtn_space, vec, name); },
      "rtn_space"_a,
      "DoF_vector"_a,
      "name"_a);

  m.def(
      "prolong",
      [](DG& coarse_dg_space, V& coarse_pressure, DG& fine_dg_space) {
        auto fine_pressure = prolong<V>(make_discrete_function(coarse_dg_space, coarse_pressure), fine_dg_space);
        return std::make_unique<V>(std::move(fine_pressure.dofs().vector()));
      },
      py::call_guard<py::gil_scoped_release>(),
      "coarse_dg_space"_a,
      "coarse_pressure"_a,
      "fine_dg_space"_a);
  m.def(
      "prolong",
      [](RTN& coarse_rtn_space, V& coarse_flux, RTN& fine_rtn_space) {
        auto fine_flux = prolong<V>(make_discrete_function(coarse_rtn_space, coarse_flux), fine_rtn_space);
        return std::make_unique<V>(std::move(fine_flux.dofs().vector()));
      },
      py::call_guard<py::gil_scoped_release>(),
      "coarse_rtn_space"_a,
      "coarse_flux"_a,
      "fine_rtn_space"_a);

  m.def(
      "assemble_SWIPDG_matrix",
      [](DG& space, XT::Functions::FunctionInterface<d>& diffusion_factor, const bool parallel) {
        const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
            XT::Common::FieldMatrix<double, 2, d>({{1., 0.}, {0., 1.}}));
        return assemble_SWIPDG_matrix(
            space, diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>(), parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "diffusion_factor"_a,
      "parallel"_a = true);
  m.def(
      "assemble_SWIPDG_matrix",
      [](DG& space,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor,
         XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
         const bool parallel) { return assemble_SWIPDG_matrix(space, diffusion_factor, diffusion_tensor, parallel); },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "diffusion_factor"_a,
      "diffusion_tensor"_a,
      "parallel"_a = true);

  m.def(
      "assemble_L2_vector",
      [](DG& space, XT::Functions::FunctionInterface<d>& force, const bool parallel) {
        return assemble_L2_vector(space, force.as_grid_function<E>(), parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "force"_a,
      "parallel"_a = true);
  m.def(
      "assemble_L2_vector",
      [](DG& space, XT::Functions::GridFunctionInterface<E>& force, const bool parallel) {
        return assemble_L2_vector(space, force, parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "force"_a,
      "parallel"_a = true);

  m.def(
      "assemble_energy_semi_product_matrix",
      [](DG& space, XT::Functions::FunctionInterface<d>& diffusion_factor, const bool parallel) {
        const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
            XT::Common::FieldMatrix<double, 2, d>({{1., 0.}, {0., 1.}}));
        return assemble_energy_semi_product_matrix(
            space, diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>(), parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "diffusion_factor"_a,
      "parallel"_a = true);
  m.def(
      "assemble_energy_semi_product_matrix",
      [](DG& space,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor,
         XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
         const bool parallel) {
        return assemble_energy_semi_product_matrix(space, diffusion_factor, diffusion_tensor, parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "diffusion_factor"_a,
      "diffusion_tensor"_a,
      "parallel"_a = true);

  m.def(
      "assemble_DG_product_matrix",
      [](DG& space, const bool parallel) {
        const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
        const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
            XT::Common::FieldMatrix<double, 2, d>({{1., 0.}, {0., 1.}}));
        return assemble_DG_product_matrix(
            space, diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>(), parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "parallel"_a = true);
  m.def(
      "assemble_DG_product_matrix",
      [](DG& space,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor,
         XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
         const bool parallel) {
        return assemble_DG_product_matrix(space, diffusion_factor, diffusion_tensor, parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "diffusion_factor"_a,
      "diffusion_tensor"_a,
      "parallel"_a = true);

  m.def(
      "compute_flux_reconstruction",
      [](GP& grid, DG& dg_space, RTN& rtn_space, XT::Functions::FunctionInterface<d>& diffusion_factor, V& dg_vec) {
        const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
            XT::Common::FieldMatrix<double, 2, d>({{1., 0.}, {0., 1.}}));
        return compute_flux_reconstruction(grid,
                                           dg_space,
                                           rtn_space,
                                           diffusion_factor.as_grid_function<E>(),
                                           diffusion_tensor.as_grid_function<E>(),
                                           dg_vec);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "dg_space"_a,
      "rtn_space"_a,
      "diffusion_factor"_a,
      "dg_DoF_vector"_a);
  m.def(
      "compute_flux_reconstruction",
      [](GP& grid,
         DG& dg_space,
         RTN& rtn_space,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor,
         XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
         V& dg_vec) {
        return compute_flux_reconstruction(grid, dg_space, rtn_space, diffusion_factor, diffusion_tensor, dg_vec);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "dg_space"_a,
      "rtn_space"_a,
      "diffusion_factor"_a,
      "diffusion_tensor"_a,
      "dg_DoF_vector"_a);

  m.def(
      "assemble_Hdiv_product_matrix",
      [](RTN& rtn_space, const bool parallel) {
        auto op = make_matrix_operator<M>(rtn_space, Stencil::element_and_intersection);
        op.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>()));
        op.append(LocalElementIntegralBilinearForm<E, d>(
            [](const auto& test_basis, const auto& ansatz_basis, const auto& /*param*/) {
              return std::max(test_basis.order() - 1, 0) + std::max(ansatz_basis.order() - 1, 0);
            },
            [](const auto& test_basis,
               const auto& ansatz_basis,
               const auto& point_in_reference_element,
               auto& result,
               const auto& /*param*/) {
              auto test_grads = test_basis.jacobians_of_set(point_in_reference_element);
              auto ansatz_grads = ansatz_basis.jacobians_of_set(point_in_reference_element);
              auto divergence = [](const auto& grad) {
                double div = 0.;
                for (size_t dd = 0; dd < d; ++dd)
                  div += grad[dd][dd];
                return div;
              };
              for (size_t ii = 0; ii < test_basis.size(); ++ii)
                for (size_t jj = 0; jj < ansatz_basis.size(); ++jj)
                  for (size_t dd = 0; dd < d; ++dd)
                    result[ii][jj] = divergence(test_grads[ii]) * divergence(ansatz_grads[jj]);
            }));
        op.assemble(parallel);
        return std::make_unique<M>(std::move(op.matrix()));
      },
      py::call_guard<py::gil_scoped_release>(),
      "dg_space"_a,
      "parallel"_a = true);

  m.def(
      "compute_local_conservation_error",
      [](GP& grid, VectorDF& flux, XT::Functions::FunctionInterface<d>& rhs, const bool parallel) {
        return compute_local_conservation_error(grid, flux, rhs.as_grid_function<E>(), parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "flux"_a,
      "rhs"_a,
      "parallel"_a = true);
  m.def(
      "compute_local_conservation_error",
      [](GP& grid, VectorDF& flux, XT::Functions::GridFunctionInterface<E>& rhs, const bool parallel) {
        return compute_local_conservation_error(grid, flux, rhs, parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "flux"_a,
      "rhs"_a,
      "parallel"_a = true);

  m.def(
      "compute_estimate",
      [](GP& grid,
         ScalarDF& pressure,
         VectorDF& flux,
         XT::Functions::FunctionInterface<d>& rhs,
         XT::Functions::FunctionInterface<d>& diffusion_factor,
         XT::Functions::FunctionInterface<d>& diffusion_factor_bar,
         XT::Functions::FunctionInterface<d>& diffusion_factor_hat,
         const double& alpha_bar,
         const double& alpha_hat,
         const double& gamma_bar,
         const int over_integrate,
         const bool& parallel) {
        const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
            XT::Common::FieldMatrix<double, 2, d>({{1., 0.}, {0., 1.}}));
        return compute_estimate(grid,
                                pressure,
                                flux,
                                rhs.as_grid_function<E>(),
                                diffusion_factor.as_grid_function<E>(),
                                diffusion_factor_bar.as_grid_function<E>(),
                                diffusion_factor_hat.as_grid_function<E>(),
                                diffusion_tensor.as_grid_function<E>(),
                                alpha_bar,
                                alpha_hat,
                                gamma_bar,
                                over_integrate,
                                parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "pressure"_a,
      "flux"_a,
      "rhs"_a,
      "diffusion_factor"_a,
      "diffusion_factor_bar"_a,
      "diffusion_factor_hat"_a,
      "alpha_bar"_a,
      "alpha_hat"_a,
      "gamma_bar"_a,
      "over_integrate"_a = 3,
      "parallel"_a = true);
  m.def(
      "compute_estimate",
      [](GP& grid,
         ScalarDF& pressure,
         VectorDF& flux,
         XT::Functions::GridFunctionInterface<E>& rhs,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor_bar,
         XT::Functions::GridFunctionInterface<E>& diffusion_factor_hat,
         XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
         const double& alpha_bar,
         const double& alpha_hat,
         const double& gamma_bar,
         const int over_integrate,
         const bool& parallel) {
        return compute_estimate(grid,
                                pressure,
                                flux,
                                rhs,
                                diffusion_factor,
                                diffusion_factor_bar,
                                diffusion_factor_hat,
                                diffusion_tensor,
                                alpha_bar,
                                alpha_hat,
                                gamma_bar,
                                over_integrate,
                                parallel);
      },
      py::call_guard<py::gil_scoped_release>(),
      "grid"_a,
      "pressure"_a,
      "flux"_a,
      "rhs"_a,
      "diffusion_factor"_a,
      "diffusion_factor_bar"_a,
      "diffusion_factor_hat"_a,
      "diffusion_tensor"_a,
      "alpha_bar"_a,
      "alpha_hat"_a,
      "gamma_bar"_a,
      "over_integrate"_a = 3,
      "parallel"_a = true);
} // PYBIND11_MODULE(...)
