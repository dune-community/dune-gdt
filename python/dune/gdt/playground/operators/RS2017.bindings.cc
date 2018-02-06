// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_ALUGRID && HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/elliptic.hh>

#include <dune/gdt/playground/operators/RS2017.hh>

using namespace Dune;
using namespace Dune::GDT::RS2017;
using XT::Grid::Layers;
using XT::Grid::Backends;
namespace py = pybind11;


            .c_str(),
        GDT::bindings::space_name<SP>::value(),
        GDT::bindings::space_name<SP>::value(),
        XT::Grid::bindings::layer_name<Layers::dd_subdomain>::value() + "_"
            + XT::Grid::bindings::backend_name<Backends::part>::value());
#ifndef NDEBUG
#endif
#ifndef NDEBUG
#endif
            class SubdomainDivergenceMatrixOperator
                : public GDT::MatrixOperatorBase<XT::LA::IstlRowMajorSparseMatrix<double>,
                                                 typename GDT::SpaceProvider<G,
                                                                             Layers::dd_subdomain,
                                                                             GDT::SpaceType::dg,
                                                                             GDT::Backends::fem,
                                                                             1,
                                                                             double,
                                                                             1>::type,
                                                 typename XT::Grid::Layer<G,
                                                                          Layers::dd_subdomain,
                                                                          Backends::part,
                                                                          XT::Grid::DD::SubdomainGrid<G>>::type,
                                                 GDT::RestrictedSpace<
                                                     typename GDT::SpaceProvider<G,
                                                                                 Layers::leaf,
                                                                                 GDT::SpaceType::rt,
                                                                                 GDT::Backends::pdelab,
                                                                                 0,
                                                                                 double,
                                                                                 G::dimension>::type,
                                                     typename XT::Grid::Layer<G,
                                                                              Layers::dd_subdomain,
                                                                              Backends::part,
                                                                              XT::Grid::DD::SubdomainGrid<G>>::type>,
                                                 double,
                                                 GDT::ChoosePattern::volume>
            {
              static_assert(XT::Grid::is_grid<G>::value, "");
              typedef GDT::MatrixOperatorBase<XT::LA::IstlRowMajorSparseMatrix<double>,
                                              typename GDT::SpaceProvider<G,
                                                                          Layers::dd_subdomain,
                                                                          GDT::SpaceType::dg,
                                                                          GDT::Backends::fem,
                                                                          1,
                                                                          double,
                                                                          1>::type,
                                              typename XT::Grid::Layer<G,
                                                                       Layers::dd_subdomain,
                                                                       Backends::part,
                                                                       XT::Grid::DD::SubdomainGrid<G>>::type,
                                              GDT::RestrictedSpace<
                                                  typename GDT::SpaceProvider<G,
                                                                              Layers::leaf,
                                                                              GDT::SpaceType::rt,
                                                                              GDT::Backends::pdelab,
                                                                              0,
                                                                              double,
                                                                              G::dimension>::type,
                                                  typename XT::Grid::Layer<G,
                                                                           Layers::dd_subdomain,
                                                                           Backends::part,
                                                                           XT::Grid::DD::SubdomainGrid<G>>::type>,
                                              double,
                                              GDT::ChoosePattern::volume>
                  BaseType;
              typedef SubdomainDivergenceMatrixOperator<G> ThisType;

            public:
              using typename BaseType::GridLayerType;
              using typename BaseType::RangeSpaceType;
              using typename BaseType::SourceSpaceType;

              typedef XT::Grid::extract_entity_t<GridLayerType> E;
              typedef XT::Grid::extract_intersection_t<GridLayerType> I;
              typedef typename G::ctype D;
              static const constexpr size_t d = G::dimension;
              typedef double R;
              typedef typename RangeSpaceType::BaseFunctionSetType RangeBasisType;

              static void bind(py::module& m)
              {
                using namespace pybind11::literals;

                py::class_<ThisType, XT::Grid::Walker<GridLayerType>> c(
                    m,
                    XT::Common::to_camel_case("RS2017_divergence_matrix_operator_subdomain_"
                                              + XT::Grid::bindings::grid_name<G>::value())
                        .c_str());
                c.def("assemble", [](ThisType& self) { self.assemble(); });
                c.def("matrix", [](ThisType& self) { return self.matrix(); });

                m.def("RS2017_make_divergence_matrix_operator_on_subdomain",
                      [](const XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
                         const ssize_t subdomain,
                         const RangeSpaceType& dg_space,
                         const SourceSpaceType& rt_space,
                         const size_t over_integrate) {
                        return new ThisType(dg_space,
                                            rt_space,
                                            dd_grid_provider.template layer<Layers::dd_subdomain, Backends::part>(
                                                XT::Common::numeric_cast<size_t>(subdomain)),
                                            over_integrate);
                      },
                      "dd_grid_provider"_a,
                      "subdomain"_a,
                      "dg_space"_a,
                      "rt_space"_a,
                      "over_integrate"_a = 2);
              } // ... bind(...)

              SubdomainDivergenceMatrixOperator(RangeSpaceType dg_space,
                                                SourceSpaceType rt_space,
                                                GridLayerType grd_lyr,
                                                const size_t over_integrate = 2)
                : BaseType(dg_space, rt_space, grd_lyr)
                , over_integrate_(over_integrate)
              {
                this->append([&](const auto& entity) {
                  const auto rt_source_basis = this->source_space().base_function_set(entity);
                  const auto dg_range_basis = this->range_space().base_function_set(entity);
                  for (size_t jj = 0; jj < rt_source_basis.size(); ++jj) {
                    const auto JJ = this->source_space().mapper().mapToGlobal(entity, jj);

                    XT::LA::CommonDenseMatrix<R> local_matrix(dg_range_basis.size(), dg_range_basis.size(), 0.);
                    XT::LA::CommonDenseVector<R> local_vector(dg_range_basis.size(), 0.);

                    typedef XT::Functions::ConstantFunction<E, D, d, R, 1> OneType;
                    const OneType one(1.);
                    const GDT::LocalVolumeIntegralOperator<GDT::LocalProductIntegrand<OneType>,
                                                           RangeBasisType,
                                                           RangeBasisType,
                                                           R>
                        local_l2_operator(one);

                    const GDT::LocalVolumeIntegralFunctional<GDT::LocalLambdaUnaryVolumeIntegrand<E, R, 1, 1>,
                                                             RangeBasisType,
                                                             R>
                    local_l2_functional([&](const auto& test_base) { return test_base.order(); },
                                        [&](const auto& test_base, const auto& xx, auto& local_vec) {
                                          const auto rt_jacs = rt_source_basis.jacobian(xx);
                                          R div = 0;
                                          for (size_t ss = 0; ss < d; ++ss)
                                            div += rt_jacs[jj][ss][ss];
                                          const auto test_vals = test_base.evaluate(xx);
                                          for (size_t ii = 0; ii < test_base.size(); ++ii)
                                            local_vec[ii] = div * test_vals[ii];
                                        });

                    local_l2_operator.apply2(dg_range_basis, dg_range_basis, local_matrix.backend());
                    local_l2_functional.apply(dg_range_basis, local_vector.backend());

                    // solve
                    XT::LA::CommonDenseVector<R> local_solution(dg_range_basis.size(), 0.);
                    try {
                      XT::LA::solve(local_matrix, local_vector, local_solution);
                    } catch (XT::LA::Exceptions::linear_solver_failed& ee) {
                      DUNE_THROW(GDT::projection_error,
                                 "Divergence projection failed because a local matrix could not be inverted!\n\n"
                                     << "This was the original error: "
                                     << ee.what());
                    }
                    for (size_t ii = 0; ii < dg_range_basis.size(); ++ii) {
                      const auto II = this->range_space().mapper().mapToGlobal(entity, ii);
                      this->matrix().set_entry(II, JJ, local_solution[ii]);
                    }
                  }
                });
              } // SubdomainDivergenceMatrixOperator(...)

              SubdomainDivergenceMatrixOperator(const ThisType&) = delete;
              SubdomainDivergenceMatrixOperator(ThisType&&) = delete;

            private:
              const size_t over_integrate_;
            }; // class SubdomainDivergenceMatrixOperator


template <class G>
            .c_str(),
        GDT::bindings::space_name<SP>::value(),
        GDT::bindings::space_name<SP>::value(),
        XT::Grid::bindings::layer_name<Layers::dd_subdomain>::value() + "_"
            + XT::Grid::bindings::backend_name<Backends::part>::value());
const auto space_name = GDT::bindings::space_name<SP>::value();
const auto grid_layer_name = XT::Grid::bindings::layer_name<Layers::dd_subdomain>::value() + "_"
                             + XT::Grid::bindings::backend_name<Backends::part>::value();
            .c_str(),
        space_name,
        space_name,
        grid_layer_name);
            GDT::bindings::internal::SystemAssembler<S, NGL>::bind(
                m,
                GDT::bindings::space_name<SP>::value(),
                GDT::bindings::space_name<SP>::value(),
                XT::Grid::bindings::layer_name<Layers::dd_subdomain_oversampled>::value() + "_"
                    + XT::Grid::bindings::backend_name<Backends::part>::value());
            PYBIND11_PLUGIN(__operators_RS2017)
            {
              using namespace pybind11::literals;

              py::module m("__operators_RS2017", "dune-gdt");
              DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.operators.RS2017");

              SwipdgPenaltySubdomainProduct<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
              SubdomainDivergenceMatrixOperator<ALU_2D_SIMPLEX_CONFORMING>::bind(m);
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
                    [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING,
                                              XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>& dd_grid_provider,
                       const ssize_t subdomain,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
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
                          min_ev = std::min(
                              min_ev,
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
                          min_ev = std::min(
                              min_ev,
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
                    [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING,
                                              XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>& dd_grid_provider,
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
                          subdomain_h =
                              std::max(subdomain_h, (subdomain_vertices[ii] - subdomain_vertices[jj]).two_norm());
                      return subdomain_h;
                    },
                    "dd_grid_provider"_a,
                    "subdomain"_a);
              m.def("RS2017_apply_l2_product",
                    [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING,
                                              XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>& dd_grid_provider,
                       const ssize_t subdomain,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& u,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& v,
                       const ssize_t over_integrate) {
                      py::gil_scoped_release DUNE_UNUSED(release);
                      return GDT::make_l2_operator(
                                 dd_grid_provider.template layer<Layers::dd_subdomain, Backends::view>(
                                     XT::Common::numeric_cast<size_t>(subdomain)),
                                 XT::Common::numeric_cast<size_t>(over_integrate))
                          ->apply2(u, v);
                    },
                    "dd_grid_provider"_a,
                    "subdomain"_a,
                    "u"_a,
                    "v"_a,
                    "over_integrate"_a = 2);

              typedef GDT::EllipticMatrixOperator<XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>,
                                                  XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>,
                                                  typename GDT::SpaceProvider<ALU_2D_SIMPLEX_CONFORMING,
                                                                              Layers::dd_subdomain,
                                                                              GDT::SpaceType::dg,
                                                                              GDT::Backends::gdt,
                                                                              1,
                                                                              double,
                                                                              1>::type,
                                                  XT::LA::IstlRowMajorSparseMatrix<double>,
                                                  typename XT::Grid::Layer<ALU_2D_SIMPLEX_CONFORMING,
                                                                           Layers::dd_subdomain,
                                                                           Backends::view>::type>
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
                    [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING,
                                              XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>& dd_grid_provider,
                       const ssize_t subdomain,
                       const typename EllipticMatrixOperatorType::SourceSpaceType& space,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda,
                       const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
                       const ssize_t over_integrate) {
                      return new EllipticMatrixOperatorType(
                          XT::Common::numeric_cast<ssize_t>(over_integrate),
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


              return m.ptr();
            }

#endif // HAVE_DUNE_ALUGRID && HAVE_DUNE_PYBINDXI
