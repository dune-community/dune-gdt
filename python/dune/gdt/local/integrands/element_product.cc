// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/local/integrands/product.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class E, size_t r = 1, class F = double>
class LocalElementProductIntegrand
{
  using G = std::decay_t<XT::Grid::extract_grid_t<E>>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalElementProductIntegrand<E, r, F>;
  using base_type = GDT::LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  template <bool scalar = (r == 1), bool anything = true>
  struct add_bind /*<scalar=true, anything>*/
  {
    static void ctor(bound_type&) {}

    static void factory(pybind11::module&, const std::string&) {}
  };

  template <bool anything>
  struct add_bind</*scalar=*/false, anything>
  {
    static void ctor(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const FieldMatrix<F, r, r>&, const std::string&>(), "weight_matrix"_a, "logging_prefix"_a = "");
      c.def(py::init<const XT::Functions::FunctionInterface<d, r, r, F>&, const std::string&>(),
            "weight_matrix_function"_a,
            "logging_prefix"_a = "",
            py::keep_alive<1, 2>());
      c.def(py::init<const XT::Functions::GridFunctionInterface<E, r, r, F>&, const std::string&>(),
            "weight_matrix_grid_function"_a,
            "logging_prefix"_a = "",
            py::keep_alive<1, 2>());
    }

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(FactoryName.c_str(),
            [](const GP&, const FieldMatrix<F, r, r>& weight_matrix, const std::string& logging_prefix) {
              return type(weight_matrix, logging_prefix);
            },
            "grid"_a,
            "weight_matrix"_a,
            "logging_prefix"_a = "");
      m.def(FactoryName.c_str(),
            [](const GP&,
               const XT::Functions::FunctionInterface<d, r, r, F>& weight_matrix_function,
               const std::string& logging_prefix) { return type(weight_matrix_function, logging_prefix); },
            "grid"_a,
            "weight_matrix_function"_a,
            "logging_prefix"_a = "",
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const XT::Functions::GridFunctionInterface<E, r, r, F>& weight_matrix_grid_function,
               const std::string& logging_prefix) { return type(weight_matrix_grid_function, logging_prefix); },
            "grid"_a,
            "weight_matrix_grid_function"_a,
            "logging_prefix"_a = "",
            py::keep_alive<0, 2>());
    }
  }; // struct add_bind<false, ...>

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_element_product_integrand",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r) + "d_bases";
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const F&, const std::string&>(), "scalar_weight"_a = F(1), "logging_prefix"_a = "");
    c.def(py::init<const XT::Functions::FunctionInterface<d, 1, 1, F>&, const std::string&>(),
          "scalar_weight_function"_a,
          "logging_prefix"_a = "",
          py::keep_alive<1, 2>());
    c.def(py::init<const XT::Functions::GridFunctionInterface<E, 1, 1, F>&, const std::string&>(),
          "scalar_weight_grid_function"_a,
          "logging_prefix"_a = "",
          py::keep_alive<1, 2>());
    add_bind<>::ctor(c);

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(FactoryName.c_str(),
          [](const GP&, const F& scalar_weight, const std::string& logging_prefix) {
            return type(scalar_weight, logging_prefix);
          },
          "grid"_a,
          "scalar_weight"_a = F(1),
          "logging_prefix"_a = "");
    m.def(FactoryName.c_str(),
          [](const GP&,
             const XT::Functions::FunctionInterface<d, 1, 1, F>& scalar_weight_function,
             const std::string& logging_prefix) { return type(scalar_weight_function, logging_prefix); },
          "grid"_a,
          "scalar_weight_function"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 2>());
    m.def(FactoryName.c_str(),
          [](const XT::Functions::GridFunctionInterface<E, 1, 1, F>& scalar_weight_grid_function,
             const std::string& logging_prefix) { return type(scalar_weight_grid_function, logging_prefix); },
          "scalar_weight_grid_function"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>());
    add_bind<>::factory(m, FactoryName);

    return c;
  } // ... bind(...)
}; // class LocalElementProductIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalElementProductIntegrand_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalElementProductIntegrand<E>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::LocalElementProductIntegrand<E, d>::bind(m);
    // add your extra dimensions here
    // ...
    LocalElementProductIntegrand_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalElementProductIntegrand_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_element_product, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_integrands_binary_element_interface");

  LocalElementProductIntegrand_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
