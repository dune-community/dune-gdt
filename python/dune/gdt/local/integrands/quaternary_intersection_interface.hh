// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_QUATERNARY_INTERSECTION_INTERFACE_HH
#define PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_QUATERNARY_INTERSECTION_INTERFACE_HH

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


template <class G,
          class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalQuaternaryIntersectionIntegrandInterface
{
protected:
  static const size_t d = G::dimension;
  using E = typename I::Entity;

public:
  using type = GDT::LocalQuaternaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using bound_type = pybind11::class_<type>;

  template <class T, typename... options>
  static void bind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "__add__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
    c.def(
        "__iadd__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
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
    std::string test_string = "";
    test_string += "_" + XT::Common::to_string(t_r) + "d";
    if (t_rC > 1)
      test_string += "x" + XT::Common::to_string(t_rC) + "d";
    if (!std::is_same<TF, double>::value)
      test_string += "_" + XT::Common::Typename<TF>::value(/*fail_wo_typeid=*/true);
    test_string += "_test_basis";
    std::string ansatz_string = "";
    ansatz_string += "_" + XT::Common::to_string(a_r) + "d";
    if (a_rC > 1)
      ansatz_string += "x" + XT::Common::to_string(a_rC) + "d";
    if (!std::is_same<AF, double>::value)
      ansatz_string += "_" + XT::Common::Typename<AF>::value(/*fail_wo_typeid=*/true);
    ansatz_string += "_ansatz_basis";
    ret += test_string;
    if (!test_string.empty() && !ansatz_string.empty())
      ret += "_x";
    ret += ansatz_string;
    ret += "_to_scalar";
    if (!std::is_same<F, double>::value)
      ret += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    return ret;
  } // ... id(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "local_quaternary_intersection_integrand")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto class_name = class_id + id(grid_id, layer_id) + "_interface";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    // static information about dims
    // ...

    bind_methods(c);

    return c;
  } // ... bind(...)
}; // class LocalQuaternaryIntersectionIntegrandInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_QUATERNARY_INTERSECTION_INTERFACE_HH
