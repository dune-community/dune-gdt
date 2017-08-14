// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <memory>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/quadraturerules.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/la/eigen-solver/eigen.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/operators/l2.hh>

//#include "RS2017.hh"

using namespace Dune;
using XT::Grid::Layers;
using XT::Grid::Backends;
namespace py = pybind11;


#if 0
template <class G, Layers layer_type, Backends layer_backend>
struct ResidualIndicatorFtimesFproduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;

  typedef GDT::RS2017::ResidualIndicator::FtimesFproduct<GL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  static std::string class_name()
  {
    return "RS2017_residual_indicator_f_times_f_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value();
  }

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "RS2017::ResidualIndicator::FtimesFproduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error&) {
    }

    m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
          [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
             const ssize_t layer_level_or_subdomain,
             const typename type::ScalarFunctionType& f,
             const ssize_t over_integrate) {
            return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                            f,
                            XT::Common::numeric_cast<size_t>(over_integrate));
          },
          "dd_grid_provider"_a,
          "layer_level"_a = -1,
          "f"_a,
          "over_integrate"_a = 2);
  } // ... bind(...)
}; // struct ResidualIndicatorFtimesFproduct


template <class G, Layers layer_type, Backends layer_backend, Layers reconstruction_layer_type>
struct ResidualIndicatorFtimesVproduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, reconstruction_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type RGL;

  typedef GDT::RS2017::ResidualIndicator::FtimesVproduct<GL, RGL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  static std::string class_name()
  {
    return "RS2017_residual_indicator_f_times_v_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + "_"
           + XT::Grid::bindings::layer_name<reconstruction_layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<Backends::view>::value();
  }

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "RS2017::ResidualIndicator::FtimesVproduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error&) {
    }

    m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
          [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
             const ssize_t layer_level_or_subdomain,
             const ssize_t reconstruction_layer_level_or_subdomain,
             const typename type::ScalarFunctionType& lambda,
             const typename type::TensorFunctionType& kappa,
             const typename type::ScalarFunctionType& f,
             const typename type::ScalarFunctionType& v,
             const ssize_t over_integrate) {
            return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                            dd_grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                                XT::Common::numeric_cast<int>(reconstruction_layer_level_or_subdomain)),
                            lambda,
                            kappa,
                            f,
                            v,
                            XT::Common::numeric_cast<size_t>(over_integrate));
          },
          "dd_grid_provider"_a,
          "layer_level"_a = -1,
          "reconstruction_layer_level"_a = -1,
          "lambda"_a,
          "kappa"_a,
          "f"_a,
          "v"_a,
          "over_integrate"_a = 2);
  } // ... bind(...)
}; // struct ResidualIndicatorFtimesVproduct


template <class G, Layers layer_type, Backends layer_backend, Layers reconstruction_layer_type>
struct ResidualIndicatorUtimesVproduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, reconstruction_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type RGL;

  typedef GDT::RS2017::ResidualIndicator::UtimesVproduct<GL, RGL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  static std::string class_name()
  {
    return "RS2017_residual_indicator_u_times_v_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + "_"
           + XT::Grid::bindings::layer_name<reconstruction_layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<Backends::view>::value();
  }

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "RS2017::ResidualIndicator::UtimesVproduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error&) {
    }

    m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
          [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
             const ssize_t layer_level_or_subdomain,
             const ssize_t reconstruction_layer_level_or_subdomain,
             const typename type::ScalarFunctionType& lambda,
             const typename type::TensorFunctionType& kappa,
             const typename type::ScalarFunctionType& u,
             const typename type::ScalarFunctionType& v,
             const ssize_t over_integrate) {
            return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                            dd_grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                                XT::Common::numeric_cast<int>(reconstruction_layer_level_or_subdomain)),
                            lambda,
                            kappa,
                            u,
                            v,
                            XT::Common::numeric_cast<size_t>(over_integrate));
          },
          "dd_grid_provider"_a,
          "layer_level"_a = -1,
          "reconstruction_layer_level"_a = -1,
          "lambda"_a,
          "kappa"_a,
          "u"_a,
          "v"_a,
          "over_integrate"_a = 2);
  } // ... bind(...)
}; // struct ResidualIndicatorUtimesVproduct
#endif // 0


PYBIND11_PLUGIN(__operators_RS2017)
{
  using namespace pybind11::literals;

  py::module m("__operators_RS2017", "dune-gdt");

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

#if HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM
//  ResidualIndicatorFtimesFproduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::part>::bind(m);
// This is not efficient: we reconstruct on the whole leaf instead of only the neighborhood, but the rt pdelab space
//                        on a dd_subdomain_oversampled grid view (which is a wrapped part) is broken, if based on
//                        a 2d simplex alugrid.
//  ResidualIndicatorFtimesVproduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::part,
//  Layers::leaf>::bind(
//      m);
//  ResidualIndicatorUtimesVproduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::part,
//  Layers::leaf>::bind(
//      m);
//  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::part, Layers::leaf>::bind(m);
#endif

  typedef typename ALU_2D_SIMPLEX_CONFORMING::template Codim<0>::Entity E;
  typedef double D;
  static const constexpr size_t d = 2;
  typedef D R;
  m.def("RS2017_residual_indicator_min_diffusion_eigenvalue",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
           const ssize_t over_int) {
          const auto over_integrate = XT::Common::numeric_cast<size_t>(over_int);
          auto subdomain_layer =
              dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
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
              XT::LA::EigenDenseMatrix<R> diffusion = local_kappa->evaluate(xx);
              diffusion *= local_lambda->evaluate(xx);
              min_ev = std::min(min_ev, XT::LA::make_eigen_solver(diffusion).min_eigenvalue());
            }
            // * and in the corners of the gigen entity.
            const auto& reference_element = ReferenceElements<D, d>::general(entity.type());
            for (int ii = 0; ii < reference_element.size(d); ++ii) {
              const auto xx = reference_element.position(ii, d);
              XT::LA::EigenDenseMatrix<R> diffusion = local_kappa->evaluate(xx);
              diffusion *= local_lambda->evaluate(xx);
              min_ev = std::min(min_ev, XT::LA::make_eigen_solver(diffusion).min_eigenvalue());
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
          auto subdomain_layer =
              dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
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
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& v,
           const ssize_t over_integrate) {
          return GDT::make_l2_operator(
                     dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
                         XT::Common::numeric_cast<size_t>(subdomain)),
                     XT::Common::numeric_cast<size_t>(over_integrate))
              ->apply2(u, v);
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "u"_a,
        "v"_a,
        "over_integrate"_a = 2);
  m.def("RS2017_diffusive_flux_indicator_apply_aa_product",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_hat,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_v,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& v,
           const ssize_t over_integrate) {
          auto subdomain_layer =
              dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
                  XT::Common::numeric_cast<size_t>(subdomain));
          XT::Grid::Walker<decltype(subdomain_layer)> walker(subdomain_layer);
          R result = 0.;
          walker.append([&](const E& entity) {
            const auto local_lambda_hat = lambda_hat.local_function(entity);
            const auto local_lambda_u = lambda_u.local_function(entity);
            const auto local_lambda_v = lambda_v.local_function(entity);
            const auto local_kappa = kappa.local_function(entity);
            const auto local_u = u.local_function(entity);
            const auto local_v = v.local_function(entity);
            const auto integrand_order = local_lambda_hat->order() + local_lambda_u->order() + local_lambda_v->order()
                                         + 3 * local_kappa->order()
                                         + size_t(std::max(ssize_t(local_u->order()) - 1, ssize_t(0)))
                                         + size_t(std::max(ssize_t(local_v->order()) - 1, ssize_t(0)));
            for (const auto& quadrature_point :
                 QuadratureRules<D, d>::rule(entity.type(), integrand_order + over_integrate)) {
              const auto xx = quadrature_point.position();
              const auto integration_factor = entity.geometry().integrationElement(xx);
              const auto quadrature_weight = quadrature_point.weight();
              XT::Common::FieldMatrix<D, d, d> diffusion_hat_inverse = local_kappa->evaluate(xx);
              XT::Common::FieldMatrix<D, d, d> diffusion_u = diffusion_hat_inverse;
              XT::Common::FieldMatrix<D, d, d> diffusion_v = diffusion_hat_inverse;
              diffusion_hat_inverse *= local_lambda_hat->evaluate(xx);
              diffusion_hat_inverse.invert(); // there is no documented way to tell if the inversion was successfull
              diffusion_u *= local_lambda_u->evaluate(xx);
              diffusion_v *= local_lambda_v->evaluate(xx);
              const auto grad_u = local_u->jacobian(xx)[0];
              const auto grad_v = local_v->jacobian(xx)[0];
              result += integration_factor * quadrature_weight
                        * ((diffusion_hat_inverse * (diffusion_u * grad_u)) * (diffusion_v * grad_v));
            }
          });
          walker.walk();
          return result;
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "lambda_hat"_a,
        "lambda_u"_a,
        "lambda_v"_a,
        "kappa"_a,
        "u"_a,
        "v"_a,
        "over_integrate"_a = 2);
  m.def("RS2017_diffusive_flux_indicator_apply_ab_product",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_hat,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d>& reconstructed_v,
           const ssize_t over_integrate) {
          auto subdomain_layer =
              dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
                  XT::Common::numeric_cast<size_t>(subdomain));
          XT::Grid::Walker<decltype(subdomain_layer)> walker(subdomain_layer);
          R result = 0.;
          walker.append([&](const E& entity) {
            const auto local_lambda_hat = lambda_hat.local_function(entity);
            const auto local_lambda_u = lambda_u.local_function(entity);
            const auto local_kappa = kappa.local_function(entity);
            const auto local_u = u.local_function(entity);
            const auto local_reconstructed_v = reconstructed_v.local_function(entity);
            const auto integrand_order = local_lambda_hat->order() + local_lambda_u->order() + 2 * local_kappa->order()
                                         + size_t(std::max(ssize_t(local_u->order()) - 1, ssize_t(0)))
                                         + local_reconstructed_v->order();
            for (const auto& quadrature_point :
                 QuadratureRules<D, d>::rule(entity.type(), integrand_order + over_integrate)) {
              const auto xx = quadrature_point.position();
              const auto integration_factor = entity.geometry().integrationElement(xx);
              const auto quadrature_weight = quadrature_point.weight();
              XT::Common::FieldMatrix<D, d, d> diffusion_hat_inverse = local_kappa->evaluate(xx);
              XT::Common::FieldMatrix<D, d, d> diffusion_u = diffusion_hat_inverse;
              diffusion_hat_inverse *= local_lambda_hat->evaluate(xx);
              diffusion_hat_inverse.invert(); // there is no documented way to tell if the inversion was successfull
              diffusion_u *= local_lambda_u->evaluate(xx);
              const auto grad_u = local_u->jacobian(xx)[0];
              const auto val_rec_v = local_reconstructed_v->evaluate(xx);
              result += integration_factor * quadrature_weight
                        * ((diffusion_hat_inverse * (diffusion_u * grad_u)) * val_rec_v);
            }
          });
          walker.walk();
          return result;
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "lambda_hat"_a,
        "lambda_u"_a,
        "kappa"_a,
        "u"_a,
        "reconstructed_v"_a,
        "over_integrate"_a = 2);
  m.def("RS2017_diffusive_flux_indicator_apply_bb_product",
        [](XT::Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING, XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>&
               dd_grid_provider,
           const ssize_t subdomain,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>& lambda_hat,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>& kappa,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d>& reconstructed_u,
           const XT::Functions::LocalizableFunctionInterface<E, D, d, R, d>& reconstructed_v,
           const ssize_t over_integrate) {
          auto subdomain_layer =
              dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
                  XT::Common::numeric_cast<size_t>(subdomain));
          XT::Grid::Walker<decltype(subdomain_layer)> walker(subdomain_layer);
          R result = 0.;
          walker.append([&](const E& entity) {
            const auto local_lambda_hat = lambda_hat.local_function(entity);
            const auto local_kappa = kappa.local_function(entity);
            const auto local_reconstructed_u = reconstructed_u.local_function(entity);
            const auto local_reconstructed_v = reconstructed_v.local_function(entity);
            const auto integrand_order = local_lambda_hat->order() + local_kappa->order()
                                         + local_reconstructed_u->order() + local_reconstructed_v->order();
            for (const auto& quadrature_point :
                 QuadratureRules<D, d>::rule(entity.type(), integrand_order + over_integrate)) {
              const auto xx = quadrature_point.position();
              const auto integration_factor = entity.geometry().integrationElement(xx);
              const auto quadrature_weight = quadrature_point.weight();
              XT::Common::FieldMatrix<D, d, d> diffusion_hat_inverse = local_kappa->evaluate(xx);
              diffusion_hat_inverse *= local_lambda_hat->evaluate(xx);
              diffusion_hat_inverse.invert(); // there is no documented way to tell if the inversion was successfull
              const auto val_rec_u = local_reconstructed_u->evaluate(xx);
              const auto val_rec_v = local_reconstructed_v->evaluate(xx);
              result += integration_factor * quadrature_weight * ((diffusion_hat_inverse * val_rec_u) * val_rec_v);
            }
          });
          walker.walk();
          return result;
        },
        "dd_grid_provider"_a,
        "subdomain"_a,
        "lambda_hat"_a,
        "kappa"_a,
        "reconstructed_u"_a,
        "reconstructed_v"_a,
        "over_integrate"_a = 2);

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = Dune::XT::Common::numeric_cast<int>(args.size());
          char** argv = Dune::XT::Common::vector_to_main_args(args);
          Dune::MPIHelper::instance(argc, argv);
#if HAVE_DUNE_FEM
          Dune::Fem::MPIManager::initialize(argc, argv);
#endif
        },
        "args"_a = std::vector<std::string>());

  m.def("_init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = std::numeric_limits<ssize_t>::max(),
        "max_debug_level"_a = std::numeric_limits<ssize_t>::max(),
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  m.def("_test_logger",
        [](const bool info, const bool debug, const bool warning) {
          auto logger = Dune::XT::Common::TimedLogger().get("dune.gdt.operators.elliptic");
          if (info)
            logger.info() << "info logging works!" << std::endl;
          if (debug)
            logger.debug() << "debug logging works!" << std::endl;
          if (warning)
            logger.warn() << "warning logging works!" << std::endl;
        },
        "info"_a = true,
        "debug"_a = true,
        "warning"_a = true);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
