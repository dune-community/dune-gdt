// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_UNARY_INTERSECTION_INTERFACE_HH
#define PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_UNARY_INTERSECTION_INTERFACE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class I, size_t r = 1, size_t rC = 1, class RF = double, class F = double>
class LocalUnaryIntersectionIntegrandInterface
{
protected:
  using G = XT::Grid::extract_grid_t<I>;
  static const size_t d = G::dimension;
  using type = GDT::LocalUnaryIntersectionIntegrandInterface<I, r, rC, RF, F>;
  using bound_type = pybind11::class_<type>;

  template <class T, typename... options>
  static void bind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def("copy_as_unary_intersection_integrand",
          [](const type& self) { return self.copy_as_unary_intersection_integrand(); });
    c.def("__add__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
    c.def("__iadd__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
    // order/evaluate
    // ...
    //    c.def("__repr__", [](const type& self) {
    //      std::stringstream ss;
    //      ss << self;
    //      return ss.str();
    //    });
  } // ... bind_methods(...)

  static std::string id(const std::string& grid_id, const std::string& layer_id)
  {
    std::string ret = grid_id;
    if (!layer_id.empty())
      ret += "_" + layer_id;
    ret += "_" + XT::Common::to_string(r) + "d";
    if (rC > 1)
      ret += "x" + XT::Common::to_string(rC) + "d";
    if (!std::is_same<RF, double>::value)
      ret += "_" + XT::Common::Typename<RF>::value(/*fail_wo_typeid=*/true);
    ret += "_test_basis";
    ret += ret;
    ret += "_to_scalar";
    if (!std::is_same<F, double>::value)
      ret += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    return ret;
  } // ... id(...)

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_unary_intersection_integrand",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto class_name = class_id + "_" + id(grid_id, layer_id) + "_interfaces";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    // static information about dims
    // ...
    bind_methods(c);

    return c;
  } // ... bind(...)
}; // class LocalUnaryIntersectionIntegrandInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_UNARY_INTERSECTION_INTERFACE_HH
