// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <memory>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include "ESV2007.hh"

using namespace Dune;
using XT::Grid::Layers;
using XT::Grid::Backends;
namespace py = pybind11;


template <class G, Layers layer_type, Backends layer_backend, Layers interpolation_layer_type = layer_type>
struct NonconformityProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, interpolation_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type IGL;

  typedef GDT::ESV2007::NonconformityProduct<GL, IGL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  template <bool is_same = (interpolation_layer_type == layer_type) && (layer_backend == Backends::part),
            bool anything = true>
  struct interpolation_layer_suffix
  {
    static std::string value()
    {
      return "";
    }
  }; // struct interpolation_layer_suffix<true, ...>

  template <bool anything>
  struct interpolation_layer_suffix<false, anything>
  {
    static std::string value()
    {
      return "_" + XT::Grid::bindings::layer_name<interpolation_layer_type>::value() + "_"
             + XT::Grid::bindings::backend_name<Backends::part>::value();
    }
  }; // struct interpolation_layer_suffix<false, ...>

  static std::string class_name()
  {
    return "ESV2007_nonconformity_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + interpolation_layer_suffix<>::value();
  }

  template <bool is_dd = (layer_type == Layers::dd_subdomain) || (layer_type == Layers::dd_subdomain_boundary)
                         || (layer_type == Layers::dd_subdomain_coupling)
                         || (layer_type == Layers::dd_subdomain_oversampled)
                         || (interpolation_layer_type == Layers::dd_subdomain)
                         || (interpolation_layer_type == Layers::dd_subdomain_boundary)
                         || (interpolation_layer_type == Layers::dd_subdomain_coupling)
                         || (interpolation_layer_type == Layers::dd_subdomain_oversampled),
            bool anything = true>
  struct factory_method
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
               const ssize_t layer_level_or_subdomain,
               const ssize_t interpolation_layer_level_or_subdomain,
               const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<IGL>>& interpolation_boundary_info,
               const typename type::ScalarFunctionType& lambda,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate) {
              return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                  XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                              dd_grid_provider.template layer<interpolation_layer_type, Backends::view>(
                                  XT::Common::numeric_cast<int>(interpolation_layer_level_or_subdomain)),
                              interpolation_boundary_info,
                              lambda,
                              kappa,
                              u,
                              v,
                              XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "dd_grid_provider"_a,
            "layer_level_or_subdomain"_a = -1,
            "interpolation_layer_level_or_subdomain"_a = -1,
            "interpolation_boundary_info"_a,
            "lambda"_a,
            "kappa"_a,
            "u"_a,
            "v"_a,
            "over_integrate"_a = 2);
    }
  }; // struct factory_method<true, ...>

  template <bool anything>
  struct factory_method<false, anything>
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider,
               const ssize_t layer_level,
               const ssize_t interpolation_layer_level,
               const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<IGL>>& interpolation_boundary_info,
               const typename type::ScalarFunctionType& lambda,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate) {
              return new type(
                  grid_provider.template layer<layer_type, layer_backend>(XT::Common::numeric_cast<int>(layer_level)),
                  grid_provider.template layer<interpolation_layer_type, Backends::view>(
                      XT::Common::numeric_cast<int>(interpolation_layer_level)),
                  interpolation_boundary_info,
                  lambda,
                  kappa,
                  u,
                  v,
                  XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "grid_provider"_a,
            "layer_level"_a = -1,
            "interpolation_layer_level"_a = -1,
            "interpolation_boundary_info"_a,
            "lambda"_a,
            "kappa"_a,
            "u"_a,
            "v"_a,
            "over_integrate"_a = 2);

      factory_method<true>::addbind(m);
    }
  }; // struct factory_method<false, ...>

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "ESV2007::NonconformityProduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
      c.def("result", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error& ee) {
    }

    factory_method<>::addbind(m);
  } // ... bind(...)
}; // struct NonconformityProduct


template <class G, Layers layer_type, Backends layer_backend, Layers reconstruction_layer_type = layer_type>
struct ResidualProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, reconstruction_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type RGL;

  typedef GDT::ESV2007::ResidualProduct<GL, RGL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  template <bool is_same = (reconstruction_layer_type == layer_type) && (layer_backend == Backends::view),
            bool anything = true>
  struct reconstruction_layer_suffix
  {
    static std::string value()
    {
      return "";
    }
  }; // struct reconstruction_layer_suffix<true, ...>

  template <bool anything>
  struct reconstruction_layer_suffix<false, anything>
  {
    static std::string value()
    {
      return "_" + XT::Grid::bindings::layer_name<reconstruction_layer_type>::value() + "_"
             + XT::Grid::bindings::backend_name<Backends::view>::value();
    }
  }; // struct reconstruction_layer_suffix<false, ...>

  static std::string class_name()
  {
    return "ESV2007_residual_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + reconstruction_layer_suffix<>::value();
  }

  template <bool is_dd = (layer_type == Layers::dd_subdomain) || (layer_type == Layers::dd_subdomain_boundary)
                         || (layer_type == Layers::dd_subdomain_coupling)
                         || (layer_type == Layers::dd_subdomain_oversampled)
                         || (reconstruction_layer_type == Layers::dd_subdomain)
                         || (reconstruction_layer_type == Layers::dd_subdomain_boundary)
                         || (reconstruction_layer_type == Layers::dd_subdomain_coupling)
                         || (reconstruction_layer_type == Layers::dd_subdomain_oversampled),
            bool anything = true>
  struct factory_method
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
               const ssize_t layer_level_or_subdomain,
               const ssize_t reconstruction_layer_level_or_subdomain,
               const typename type::ScalarFunctionType& lambda,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& f,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate,
               const double& poincare_constant) {
              return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                  XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                              dd_grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                                  XT::Common::numeric_cast<int>(reconstruction_layer_level_or_subdomain)),
                              lambda,
                              kappa,
                              f,
                              u,
                              v,
                              poincare_constant,
                              XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "dd_grid_provider"_a,
            "layer_level"_a = -1,
            "reconstruction_layer_level"_a = -1,
            "lambda"_a,
            "kappa"_a,
            "f"_a,
            "u"_a,
            "v"_a,
            "over_integrate"_a = 2,
            "poincare_constant"_a = 1.0 / (M_PIl * M_PIl));
    }
  }; // struct factory_method<true, ...>

  template <bool anything>
  struct factory_method<false, anything>
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider,
               const ssize_t layer_level,
               const ssize_t reconstruction_layer_level,
               const typename type::ScalarFunctionType& lambda,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& f,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate,
               const double& poincare_constant) {
              return new type(
                  grid_provider.template layer<layer_type, layer_backend>(XT::Common::numeric_cast<int>(layer_level)),
                  grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                      XT::Common::numeric_cast<int>(reconstruction_layer_level)),
                  lambda,
                  kappa,
                  f,
                  u,
                  v,
                  poincare_constant,
                  XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "grid_provider"_a,
            "layer_level"_a = -1,
            "reconstruction_layer_level"_a = -1,
            "lambda"_a,
            "kappa"_a,
            "f"_a,
            "u"_a,
            "v"_a,
            "over_integrate"_a = 2,
            "poincare_constant"_a = 1.0 / (M_PIl * M_PIl));

      factory_method<true>::addbind(m);
    }
  }; // struct factory_method<false, ...>

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "ESV2007::ResidualProduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
      c.def("result", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error& ee) {
    }

    factory_method<>::addbind(m);
  } // ... bind(...)
}; // struct ResidualProduct


template <class G, Layers layer_type, Backends layer_backend, Layers reconstruction_layer_type = layer_type>
struct DiffusiveFluxProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, reconstruction_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type RGL;

  typedef GDT::ESV2007::DiffusiveFluxProduct<GL, RGL> type;
  typedef py::class_<type, XT::Grid::Walker<GL>> bound_type;

  template <bool is_same = (reconstruction_layer_type == layer_type) && (layer_backend == Backends::view),
            bool anything = true>
  struct reconstruction_layer_suffix
  {
    static std::string value()
    {
      return "";
    }
  }; // struct reconstruction_layer_suffix<true, ...>

  template <bool anything>
  struct reconstruction_layer_suffix<false, anything>
  {
    static std::string value()
    {
      return "_" + XT::Grid::bindings::layer_name<reconstruction_layer_type>::value() + "_"
             + XT::Grid::bindings::backend_name<Backends::view>::value();
    }
  }; // struct reconstruction_layer_suffix<false, ...>

  static std::string class_name()
  {
    return "ESV2007_diffusive_flux_product";
  }

  static std::string layer_suffix()
  {
    return XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + reconstruction_layer_suffix<>::value();
  }

  template <bool is_dd = (layer_type == Layers::dd_subdomain) || (layer_type == Layers::dd_subdomain_boundary)
                         || (layer_type == Layers::dd_subdomain_coupling)
                         || (layer_type == Layers::dd_subdomain_oversampled)
                         || (reconstruction_layer_type == Layers::dd_subdomain)
                         || (reconstruction_layer_type == Layers::dd_subdomain_boundary)
                         || (reconstruction_layer_type == Layers::dd_subdomain_coupling)
                         || (reconstruction_layer_type == Layers::dd_subdomain_oversampled),
            bool anything = true>
  struct factory_method
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

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
    }
  }; // struct factory_method<true, ...>

  template <bool anything>
  struct factory_method<false, anything>
  {
    static void addbind(py::module& m)
    {
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name() + "_" + layer_suffix()).c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider,
               const ssize_t layer_level,
               const ssize_t reconstruction_layer_level,
               const typename type::ScalarFunctionType& lambda,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate) {
              return new type(
                  grid_provider.template layer<layer_type, layer_backend>(XT::Common::numeric_cast<int>(layer_level)),
                  grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                      XT::Common::numeric_cast<int>(reconstruction_layer_level)),
                  lambda,
                  kappa,
                  u,
                  v,
                  XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "grid_provider"_a,
            "layer_level"_a = -1,
            "reconstruction_layer_level"_a = -1,
            "lambda"_a,
            "kappa"_a,
            "u"_a,
            "v"_a,
            "over_integrate"_a = 2);

      factory_method<true>::addbind(m);
    }
  }; // struct factory_method<false, ...>

  static void bind(py::module& m)
  {
    using namespace pybind11::literals;

    try { // we might not be the first ones to add this type
      bound_type c(m,
                   XT::Common::to_camel_case(class_name() + "_" + XT::Grid::bindings::grid_name<G>::value() + "_"
                                             + layer_suffix())
                       .c_str(),
                   "ESV2007::DiffusiveFluxProduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
      c.def("result", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error& ee) {
    }

    factory_method<>::addbind(m);
  } // ... bind(...)
}; // struct DiffusiveFluxProduct


PYBIND11_PLUGIN(__operators_ESV2007)
{
  using namespace pybind11::literals;

  py::module m("__operators_ESV2007", "dune-gdt");

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

#if HAVE_DUNE_ALUGRID
  NonconformityProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  NonconformityProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view>::bind(m);
  NonconformityProduct<ALU_2D_SIMPLEX_CONFORMING,
                       Layers::dd_subdomain,
                       Backends::view,
                       Layers::dd_subdomain_oversampled>::bind(m);
  ResidualProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  // This is not efficient: we reconstruct on the whole leaf instead of only the neighborhood, but the rt space
  //                        on a dd_subdomain_oversampled grid view is broken, if based on
  //                        a 2d simplex alugrid.
  ResidualProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view, Layers::leaf>::bind(m);
  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  // s.a.
  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view, Layers::leaf>::bind(m);
#endif

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
