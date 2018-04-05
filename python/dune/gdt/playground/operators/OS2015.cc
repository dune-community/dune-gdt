// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#include <memory>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <python/dune/gdt/playground/operators/OS2015.hh>

using namespace Dune;
using XT::Grid::Layers;
using XT::Grid::Backends;
namespace py = pybind11;


template <class G, Layers layer_type, Backends layer_backend, Layers reconstruction_layer_type = layer_type>
struct ResidualProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef
      typename XT::Grid::Layer<G, reconstruction_layer_type, Backends::view, XT::Grid::DD::SubdomainGrid<G>>::type RGL;

  typedef GDT::OS2015::ResidualProduct<GL, RGL> type;
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
    return "OS2015_residual_product";
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
               const typename type::ScalarFunctionType& lambda_hat,
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
                              lambda_hat,
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
            "lambda_hat"_a,
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
               const typename type::ScalarFunctionType& lambda_hat,
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
                  lambda_hat,
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
            "lambda_hat"_a,
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
                   "OS2015::ResidualProduct");
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

  typedef GDT::OS2015::DiffusiveFluxProduct<GL, RGL> type;
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
    return "OS2015_diffusive_flux_product";
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
               const typename type::ScalarFunctionType& lambda_hat,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate) {
              return new type(dd_grid_provider.template layer<layer_type, layer_backend>(
                                  XT::Common::numeric_cast<int>(layer_level_or_subdomain)),
                              dd_grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                                  XT::Common::numeric_cast<int>(reconstruction_layer_level_or_subdomain)),
                              lambda,
                              lambda_hat,
                              kappa,
                              u,
                              v,
                              XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "dd_grid_provider"_a,
            "layer_level"_a = -1,
            "reconstruction_layer_level"_a = -1,
            "lambda"_a,
            "lambda_hat"_a,
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
               const typename type::ScalarFunctionType& lambda_hat,
               const typename type::TensorFunctionType& kappa,
               const typename type::ScalarFunctionType& u,
               const typename type::ScalarFunctionType& v,
               const ssize_t over_integrate) {
              return new type(
                  grid_provider.template layer<layer_type, layer_backend>(XT::Common::numeric_cast<int>(layer_level)),
                  grid_provider.template layer<reconstruction_layer_type, Backends::view>(
                      XT::Common::numeric_cast<int>(reconstruction_layer_level)),
                  lambda,
                  lambda_hat,
                  kappa,
                  u,
                  v,
                  XT::Common::numeric_cast<size_t>(over_integrate));
            },
            "grid_provider"_a,
            "layer_level"_a = -1,
            "reconstruction_layer_level"_a = -1,
            "lambda"_a,
            "lambda_hat"_a,
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
                   "OS2015::DiffusiveFluxProduct");
      c.def("apply2", [](type& self) { return self.apply2(); });
      c.def("result", [](type& self) { return self.apply2(); });
    } catch (std::runtime_error& ee) {
    }

    factory_method<>::addbind(m);
  } // ... bind(...)
}; // struct DiffusiveFluxProduct


PYBIND11_MODULE(__operators_OS2015, m)
{
  using namespace pybind11::literals;


#if HAVE_DUNE_ALUGRID
  // This is not efficient: we reconstruct on the whole leaf instead of only the neighborhood, but the rt space
  //                        on a dd_subdomain_oversampled grid view (which is a wrapped part) is broken, if based on
  //                        a 2d simplex alugrid.
  ResidualProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view, Layers::leaf>::bind(m);
  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  // s.a.
  DiffusiveFluxProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view, Layers::leaf>::bind(m);
#endif
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.operators.elliptic");
}
