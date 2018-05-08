// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_OPERATORS_L2_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_L2_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.hh>
#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include <python/dune/gdt/operators/base.hh>
#include <dune/gdt/operators/l2.hh>

namespace Dune {
namespace GDT {
namespace bindings {
namespace internal {


template <class R, class M>
class L2MatrixOperator
{
public:
  typedef GDT::L2MatrixOperator<R, M> type;
  using bound_type = pybind11::class_<type, typename bindings::MatrixOperatorBase<type>::BaseType>;

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& space_name,
                         const std::string& container_name,
                         const std::string& grid_layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const std::string class_name = "l2_matrix_operator";
    const auto ClassName = XT::Common::to_camel_case(class_name + "_" + space_name + "_" + container_name);

    XT::Common::bindings::try_register(m, [&](pybind11::module& mod) {
      MatrixOperatorBase<type>::bind(mod, ClassName, space_name, space_name, grid_layer_name);
    });

    bound_type c(m, ClassName.c_str());
    c.def("assemble", [](type& self) { self.assemble(); });
    c.def("matrix", [](type& self) { return self.matrix(); });

    m.def(std::string("make_" + class_name + "_" + container_name).c_str(),
          [](const R& space, const size_t over_integrate) {
            return make_l2_matrix_operator<M>(space, over_integrate).release();
          },
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>());

    m.def(std::string("make_" + class_name).c_str(),
          [](M& matrix, const R& space, const size_t over_integrate) {
            return make_l2_matrix_operator(matrix, space, over_integrate).release();
          },
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class L2MatrixOperator


} // namespace internal


template <class G,
          XT::Grid::Layers layer_type,
          GDT::SpaceType space_type,
          GDT::Backends space_backend,
          int p,
          size_t r,
          XT::LA::Backends la_backend>
class L2MatrixOperator
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef GDT::SpaceProvider<G, layer_type, space_type, space_backend, p, double, r, 1> RP;
  typedef typename RP::type R;
  typedef typename XT::LA::Container<double, la_backend>::MatrixType M;

  typedef internal::L2MatrixOperator<R, M> binder;

public:
  typedef typename binder::type type;
  typedef typename binder::bound_type bound_type;

public:
  static bound_type bind(pybind11::module& m)
  {
    const auto grid_layer_name = XT::Grid::bindings::layer_name<layer_type>::value()
                                 + XT::Grid::bindings::backend_name<XT::Grid::Backends::view>::value();
    return binder::bind(m, space_name<RP>::value(), XT::LA::bindings::container_name<M>::value(), grid_layer_name);
  }
}; // class L2MatrixOperator


template <class G,
          XT::Grid::Layers layer_type,
          XT::Grid::Backends layer_backend,
          size_t range_r = 1,
          size_t range_rC = 1,
          size_t source_r = range_r,
          size_t source_rC = range_rC>
class L2LocalizableProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");

  template <bool is_dd = layer_type == XT::Grid::Layers::dd_subdomain
                         || layer_type == XT::Grid::Layers::dd_subdomain_boundary
                         || layer_type == XT::Grid::Layers::dd_subdomain_coupling
                         || layer_type == XT::Grid::Layers::dd_subdomain_oversampled,
            bool anything = true>
  struct GridLayer
  {
    typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type type;
  };

  template <bool anything>
  struct GridLayer<false, anything>
  {
    typedef typename XT::Grid::Layer<G, layer_type, layer_backend>::type type;
  };

  typedef typename GridLayer<>::type GL;
  typedef XT::Grid::extract_entity_t<GL> E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, double, range_r, range_rC> R;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, double, source_r, source_rC> S;

public:
  typedef GDT::L2LocalizableProduct<GL, R, S> type;
  typedef pybind11::class_<type, XT::Grid::Walker<GL>> bound_type;

private:
  static std::string class_name()
  {
    return "l2_localizable_product_on_" + XT::Grid::bindings::layer_name<layer_type>::value() + "_"
           + XT::Grid::bindings::backend_name<layer_backend>::value() + "_for_" + XT::Common::to_string(range_r) + "x"
           + XT::Common::to_string(range_rC) + "_range_times_" + XT::Common::to_string(source_r) + "x"
           + XT::Common::to_string(source_rC) + "_source";
  }

  template <bool is_dd = layer_type == XT::Grid::Layers::dd_subdomain
                         || layer_type == XT::Grid::Layers::dd_subdomain_boundary
                         || layer_type == XT::Grid::Layers::dd_subdomain_coupling
                         || layer_type == XT::Grid::Layers::dd_subdomain_oversampled,
            bool anything = true>
  struct FactoryMethods
  {
    static void addbind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name()).c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider,
               const int level_or_subdomain,
               const R& range,
               const S& source,
               const size_t over_integrate) {
              return make_l2_localizable_product(
                         grid_provider.template layer<layer_type, layer_backend>(level_or_subdomain),
                         range,
                         source,
                         over_integrate)
                  .release();
            },
            "grid_provider"_a,
            "level_or_subdomain"_a,
            "range"_a,
            "source"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 3>(),
            py::keep_alive<0, 4>());
    }
  }; // struct FactoryMethods<true, ...>

  template <bool anything>
  struct FactoryMethods<false, anything>
  {
    static void addbind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(std::string("make_" + class_name()).c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider,
               const int level,
               const R& range,
               const S& source,
               const size_t over_integrate) {
              return make_l2_localizable_product(
                         grid_provider.template layer<layer_type, layer_backend>(level), range, source, over_integrate)
                  .release();
            },
            "grid_provider"_a,
            "level"_a,
            "range"_a,
            "source"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 3>(),
            py::keep_alive<0, 4>());

      FactoryMethods<true>::addbind(m);
    }
  }; // struct FactoryMethods<false, ...>

public:
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    FactoryMethods<>::addbind(m); // needs to come first, as the code below may fail

    const std::string cn = class_name();
    bound_type c(m, XT::Common::to_camel_case(cn).c_str());
    c.def("apply2", [](type& self) { return self.apply2(); });
    c.def("result", [](type& self) { return self.result(); });

    return c;
  } // ... bind(...)
}; // class L2LocalizableProduct


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_OPERATORS_L2_BINDINGS_HH
