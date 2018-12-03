// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#include "config.h"

#if HAVE_DUNE_ALUGRID && HAVE_DUNE_PYBINDXI

#  include <dune/common/parallel/mpihelper.hh>

#  include <dune/pybindxi/pybind11.h>
#  include <dune/pybindxi/stl.h>

#  include <python/dune/xt/common/bindings.hh>
#  include <python/dune/gdt/shared.hh>

#  include <dune/xt/la/eigen-solver.hh>

#  include <dune/gdt/operators/elliptic.hh>

#  include <python/dune/gdt/playground/operators/RS2017.hh>

using namespace Dune;
using namespace Dune::GDT::RS2017;
using XT::Grid::Backends;
using XT::Grid::Layers;
namespace py = pybind11;


PYBIND11_MODULE(__operators_RS2017, m)
{
  using namespace pybind11::literals;

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

  SwipdgPenaltySubdomainProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
  HdivSemiProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
  DiffusiveFluxAaProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
  DiffusiveFluxAbProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
  DiffusiveFluxBbProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
  ResidualPartFunctional<ALU_2D_SIMPLEX_CONFORMING>::bind(m);

  bind_neighborhood_reconstruction<ALU_2D_SIMPLEX_CONFORMING>(m);
  bind_neighborhood_discretization<ALU_2D_SIMPLEX_CONFORMING>(m);

  typedef typename ALU_2D_SIMPLEX_CONFORMING::template Codim<0>::Entity E;
  typedef double D;
  static const constexpr size_t d = 2;
  typedef double R;

  typedef XT::LA::IstlDenseVector<R> V;
  typedef XT::LA::IstlRowMajorSparseMatrix<R> M;

  m.def("RS2017_residual_indicator_min_diffusion_eigenvalue",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::GridFunctionInterface<E, D, d, R, 1>& lambda,
           const XT::Functions::GridFunctionInterface<E, D, d, R, d, d>& kappa,
           const ssize_t over_int) {
          py::gil_scoped_release DUNE_UNUSED(release);
          const auto over_integrate = XT::Common::numeric_cast<size_t>(over_int);
          auto subdomain_layer = dd_grid_provider.template layer<Layers::dd_subdomain, Backends::view>(
              XT::Common::numeric_cast<size_t>(subdomain));
          typedef decltype(subdomain_layer) GL;
          XT::Grid::Walker<GL> walker(subdomain_layer);
          double min_ev = std::numeric_limits<double>::max();
          walker.append([&](const E& entity) {
            const auto local_lambda = lambda.local_function(entity);
            const auto local_kappa = kappa.local_function(entity);
            // To find the minimum of a function we evaluate it
            // * in all quadrature points of a quadrature which would integrate such a function exactly
            for (const auto& quadrature_point : QuadratureRules<D, d>::rule(
                     entity.type(), local_lambda->order() + local_kappa->order() + over_integrate)) {
              const auto xx = quadrature_point.position();
              auto diffusion = local_kappa->evaluate(xx);
              diffusion *= local_lambda->evaluate(xx);
              min_ev = std::min(min_ev,
                                XT::LA::make_eigen_solver(
                                    diffusion, XT::Common::Configuration{{"assert_positive_eigenvalues"}, {1e-15}})
                                    .min_eigenvalues(1)
                                    .at(0));
            }
            // * and in the corners of the gigen entity.
            const auto& reference_element = ReferenceElements<D, d>::general(entity.type());
            for (int ii = 0; ii < reference_element.size(d); ++ii) {
              const auto xx = reference_element.position(ii, d);
              auto diffusion = local_kappa->evaluate(xx);
              diffusion *= local_lambda->evaluate(xx);
              min_ev = std::min(min_ev,
                                XT::LA::make_eigen_solver(
                                    diffusion, XT::Common::Configuration{{"assert_positive_eigenvalues"}, {1e-15}})
                                    .min_eigenvalues(1)
                                    .at(0));
            }
          });
          walker.walk();
          return min_ev;
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "lambda"_a,
        "kappa"_a,
        "over_integrate"_a = 2);
  m.def("RS2017_residual_indicator_subdomain_diameter",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain) {
          py::gil_scoped_release DUNE_UNUSED(release);
          auto subdomain_layer = dd_grid_provider.template layer<Layers::dd_subdomain, Backends::view>(
              XT::Common::numeric_cast<size_t>(subdomain));
          typedef decltype(subdomain_layer) GL;
          XT::Grid::Walker<GL> walker(subdomain_layer);
          std::vector<FieldVector<D, d>> subdomain_vertices;
          walker.append([&](const E& entity) {
            for (size_t cc = 0; cc < entity.subEntities(d); ++cc)
              subdomain_vertices.emplace_back(entity.template subEntity<d>(cc).geometry().center());
          });
          walker.walk();
          R subdomain_h = std::numeric_limits<R>::min();
          for (size_t ii = 0; ii < subdomain_vertices.size(); ++ii)
            for (size_t jj = ii + 1; jj < subdomain_vertices.size(); ++jj)
              subdomain_h = std::max(subdomain_h, (subdomain_vertices[ii] - subdomain_vertices[jj]).two_norm());
          return subdomain_h;
        },
        "dd_grid_provider"_a,
        "subdomain"_a);
  m.def("RS2017_apply_l2_product",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::GridFunctionInterface<E, D, d, R, 1>& u,
           const XT::Functions::GridFunctionInterface<E, D, d, R, 1>& v,
           const ssize_t over_integrate) {
          py::gil_scoped_release DUNE_UNUSED(release);
          return GDT::make_l2_operator(dd_grid_provider.template layer<Layers::dd_subdomain, Backends::view>(
                                           XT::Common::numeric_cast<size_t>(subdomain)),
                                       XT::Common::numeric_cast<size_t>(over_integrate))
              ->apply2(u, v);
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "u"_a,
        "v"_a,
        "over_integrate"_a = 2);

  typedef GDT::EllipticMatrixOperator<
      XT::Functions::GridFunctionInterface<E, D, d, R, 1>,
      XT::Functions::GridFunctionInterface<E, D, d, R, d, d>,
      typename GDT::SpaceProvider<ALU_2D_SIMPLEX_CONFORMING,
                                  Layers::dd_subdomain,
                                  GDT::SpaceType::dg,
                                  GDT::Backends::gdt,
                                  1,
                                  double,
                                  1>::type,
      XT::LA::IstlRowMajorSparseMatrix<double>,
      typename XT::Grid::Layer<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view>::type>
      EllipticMatrixOperatorType;
  try { // we might not be the first to add this
    py::class_<EllipticMatrixOperatorType,
               GDT::SystemAssembler<typename EllipticMatrixOperatorType::SourceSpaceType,
                                    typename EllipticMatrixOperatorType::GridLayerType>>
        elliptic_matrix_operator(m, "EllipticMatrixOperatorNeighborhood");
    elliptic_matrix_operator.def("matrix", [](EllipticMatrixOperatorType& self) { return self.matrix(); });
  } catch (std::runtime_error&) {
  }
  m.def("RS2017_make_elliptic_matrix_operator_on_subdomain",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const typename EllipticMatrixOperatorType::SourceSpaceType& space,
           const XT::Functions::GridFunctionInterface<E, D, d, R, 1>& lambda,
           const XT::Functions::GridFunctionInterface<E, D, d, R, d, d>& kappa,
           const ssize_t over_integrate) {
          return new EllipticMatrixOperatorType(XT::Common::numeric_cast<ssize_t>(over_integrate),
                                                lambda,
                                                kappa,
                                                space,
                                                dd_grid_provider.template layer<Layers::dd_subdomain, Backends::view>(
                                                    XT::Common::numeric_cast<int>(subdomain)));
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "space"_a,
        "lambda"_a,
        "kappa"_a,
        "over_integrate"_a = 2);

  add_initialization(m, "dune.gdt.operators.elliptic");
}

#endif // HAVE_DUNE_ALUGRID && HAVE_DUNE_PYBINDXI
