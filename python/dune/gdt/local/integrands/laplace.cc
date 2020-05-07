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

#include <dune/gdt/local/integrands/laplace.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class E, size_t r = 1, class F = double>
class LocalLaplaceIntegrand
{
  using G = std::decay_t<XT::Grid::extract_grid_t<E>>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalLaplaceIntegrand<E, r, F>;
  using base_type = GDT::LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  template <bool oned = (d == 1), bool anything = true>
  struct add_bind /*<true, anything>*/
  {
    static void ctor(bound_type&) {}

    static void factory(pybind11::module&, const std::string&) {}
  };

  template <bool anything>
  struct add_bind<false, anything>
  {
    static void ctor(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const FieldMatrix<F, d, d>&>(), "constant_diffusion_tensor"_a);
      c.def(py::init<const XT::Functions::FunctionInterface<d, d, d, F>&>(),
            "diffusion_tensor_function"_a,
            py::keep_alive<1, 2>());
    }

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(FactoryName.c_str(),
            [](const GP&, const FieldMatrix<F, d, d>& constant_diffusion_tensor) {
              return type(constant_diffusion_tensor);
            },
            "unused_grid_to_select_type"_a,
            "constant_diffusion_tensor"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const XT::Functions::FunctionInterface<d, d, d, F>& diffusion_tensor_function) {
              return type(diffusion_tensor_function);
            },
            "unused_grid_to_select_type"_a,
            "diffusion_tensor_function"_a,
            py::keep_alive<0, 2>());
    }
  }; // struct add_bind<false, ...>

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_laplace_integrand",
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
    c.def(py::init<const F&>(), "constant_scalar_diffusion"_a);
    c.def(py::init<const XT::Functions::FunctionInterface<d, 1, 1, F>&>(),
          "scalar_diffusion_function"_a,
          py::keep_alive<1, 2>());
    add_bind<>::ctor(c);
    c.def(py::init<const XT::Functions::GridFunctionInterface<E, 1, 1, F>&>(),
          "scalar_diffusion_grid_function"_a,
          py::keep_alive<1, 2>());
    c.def(py::init<const XT::Functions::GridFunctionInterface<E, d, d, F>&>(),
          "diffusion_tensor_grid_function"_a,
          py::keep_alive<1, 2>());

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(FactoryName.c_str(),
          [](const GP&, const F& constant_scalar_diffusion) { return type(constant_scalar_diffusion); },
          "unused_grid_to_select_type"_a,
          "constant_scalar_diffusion"_a);
    add_bind<>::factory(m, FactoryName);
    m.def(FactoryName.c_str(),
          [](const GP&, const XT::Functions::FunctionInterface<d, 1, 1, F>& scalar_diffusion_function) {
            return type(scalar_diffusion_function);
          },
          "unused_grid_to_select_type"_a,
          "scalar_diffusion_function"_a,
          py::keep_alive<0, 2>());
    m.def(FactoryName.c_str(),
          [](const XT::Functions::GridFunctionInterface<E, 1, 1, F>& scalar_diffusion_grid_function) {
            return type(scalar_diffusion_grid_function);
          },
          "scalar_diffusion_grid_function"_a,
          py::keep_alive<0, 1>());
    m.def(FactoryName.c_str(),
          [](const XT::Functions::GridFunctionInterface<E, d, d, F>& diffusion_tensor_grid_function) {
            return type(diffusion_tensor_grid_function);
          },
          "diffusion_tensor_grid_function"_a,
          py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)
}; // class LocalLaplaceIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalLaplaceIntegrand_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalLaplaceIntegrand<E>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::LocalLaplaceIntegrand<E, d>::bind(m);
    // add your extra dimensions here
    // ...
    LocalLaplaceIntegrand_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalLaplaceIntegrand_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_laplace, m)
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

  LocalLaplaceIntegrand_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
