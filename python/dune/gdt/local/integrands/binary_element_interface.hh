// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_ELEMENT_INTERFACE_HH
#define PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_ELEMENT_INTERFACE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "conversion.hh"


namespace Dune {
namespace GDT {
namespace bindings {


template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryElementIntegrandInterface
{
protected:
  using G = XT::Grid::extract_grid_t<E>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using bound_type = pybind11::class_<type>;

private:
  template <bool square_matrix = (a_rC == a_r && a_r > 1), bool scalar = (a_rC == 1 && a_r == 1), bool anything = true>
  struct with_ansatz /*the general vector or non-square matrix case*/
  {
    static void addbind(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def("with_ansatz",
            [](type& self, const typename XT::Functions::RangeTypeSelector<F, a_r, a_rC>::type& constant_value) {
              return self.with_ansatz(constant_value);
            },
            "constant_value"_a,
            py::keep_alive<0, 1>());
      c.def("with_ansatz",
            [](type& self, const XT::Functions::FunctionInterface<d, a_r, a_rC, F>& function) {
              return self.with_ansatz(function);
            },
            "scalar_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
      c.def("with_ansatz",
            [](type& self, const XT::Functions::GridFunctionInterface<E, a_r, a_rC, F>& grid_function) {
              return self.with_ansatz(grid_function);
            },
            "scalar_grid_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    }
  }; // struct with_ansatz<true, false, ...>

  template <bool anything>
  struct with_ansatz</*square_matrix=*/true, false, anything>
  {
    static void addbind(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // we have the scalar case ...
      with_ansatz<false, true, anything>::addbind(c);

      // .. and the tensor case
      c.def("with_ansatz",
            [](type& self, const FieldMatrix<F, a_r, a_rC>& square_matrix) { return self.with_ansatz(square_matrix); },
            "square_matrix"_a,
            py::keep_alive<0, 1>());
      c.def("with_ansatz",
            [](type& self, const XT::Functions::FunctionInterface<d, a_r, a_rC, F>& square_matrix_valued_function) {
              return self.with_ansatz(square_matrix_valued_function);
            },
            "square_matrix_valued_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
      c.def("with_ansatz",
            [](type& self,
               const XT::Functions::GridFunctionInterface<E, a_r, a_rC, F>& square_matrix_valued_grid_function) {
              return self.with_ansatz(square_matrix_valued_grid_function);
            },
            "square_matrix_valued_grid_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    }
  }; // struct with_ansatz<true, false, ...>

  template <bool anything>
  struct with_ansatz<false, /*scalar=*/true, anything>
  {
    static void addbind(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def("with_ansatz",
            [](type& self, const F& constant_scalar) { return self.with_ansatz(constant_scalar); },
            "constant_scalar"_a,
            py::keep_alive<0, 1>());
      c.def("with_ansatz",
            [](type& self, const XT::Functions::FunctionInterface<d, 1, 1, F>& scalar_function) {
              return self.with_ansatz(scalar_function);
            },
            "scalar_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
      c.def("with_ansatz",
            [](type& self, const XT::Functions::GridFunctionInterface<E, 1, 1, F>& scalar_grid_function) {
              return self.with_ansatz(scalar_grid_function);
            },
            "scalar_grid_function"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    }
  }; // struct with_ansatz<true, false, ...>

protected:
  template <class T, typename... options>
  static void bind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def("copy_as_binary_element_integrand", [](const type& self) { return self.copy_as_binary_element_integrand(); });
    c.def("__add__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
    c.def("__iadd__", [](type& self, type& other) { return self + other; }, "other"_a, py::is_operator());
    // order/evaluate
    // ...
    //    c.def("__repr__", [](const type& self) {
    //      std::stringstream ss;
    //      ss << self;
    //      return ss.str();
    //    });

    // conversion to unary
    with_ansatz<>::addbind(c);
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

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_binary_element_integrand",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
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
}; // class LocalBinaryElementIntegrandInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_ELEMENT_INTERFACE_HH
