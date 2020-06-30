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
#include <dune/gdt/local/functionals/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class E, size_t r = 1, size_t rC = 1, class RF = double, class F = double>
class LocalElementFunctionalInterface
{
  using G = XT::Grid::extract_grid_t<E>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalElementFunctionalInterface<E, r, rC, RF, F>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_element_functional",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    std::string test_string = "";
    test_string += "_" + XT::Common::to_string(r) + "d";
    if (rC > 1)
      test_string += "x" + XT::Common::to_string(rC) + "d";
    if (!std::is_same<RF, double>::value)
      test_string += "_" + XT::Common::Typename<RF>::value(/*fail_wo_typeid=*/true);
    test_string += "_test_basis";
    class_name += test_string;
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    class_name += "_interface";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    // static information about dims
    // ...
    c.def("copy", [](const type& self) { return self.copy(); });
    // apply2
    // ...
    //    c.def("__repr__", [](const type& self) {
    //      std::stringstream ss;
    //      ss << self;
    //      return ss.str();
    //    });
    return c;
  } // ... bind(...)
}; // class LocalElementFunctionalInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalElementFunctionalInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;
  using F = double;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalElementFunctionalInterface<E>::bind(m);
    if (d > 1) {
      Dune::GDT::bindings::LocalElementFunctionalInterface<E, d, 1, F, F>::bind(m);
      Dune::GDT::bindings::LocalElementFunctionalInterface<E, d, d, F, F>::bind(m);
    }
    // add your extra dimensions here
    // ...
    LocalElementFunctionalInterface_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalElementFunctionalInterface_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_functionals_element_interface, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  LocalElementFunctionalInterface_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
