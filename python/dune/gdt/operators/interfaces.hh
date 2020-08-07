// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef PYTHON_DUNE_GDT_OPERATORS_INTERFACES_HH
#define PYTHON_DUNE_GDT_OPERATORS_INTERFACES_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class GV, size_t s_r = 1, size_t r_r = s_r>
class OperatorInterface
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  using type = GDT::OperatorInterface<M, GV, s_r, 1, r_r, 1>;
  using V = typename type::VectorType;
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using F = typename type::FieldType;
  using Mop = typename type::MatrixOperatorType;
  using CLop = typename type::ConstLincombOperatorType;
  using Lop = typename type::LincombOperatorType;

public:
  using bound_type = pybind11::class_<type>;

  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "assemble", [](T& self, const bool parallel) { self.assemble(parallel); }, "parallel"_a = false);
    c.def(
        "apply",
        [](T& self, const V& source, V& range, const XT::Common::Parameter& param) {
          self.apply(source, range, param);
        },
        "source"_a,
        "range"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply",
        [](T& self, const V& source, const XT::Common::Parameter& param) { return self.apply(source, param); },
        "source"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply2",
        [](T& self, const V& range, const V& source, const XT::Common::Parameter& param) {
          self.apply2(range, source, param);
        },
        "range"_a,
        "source"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def("invert_options", [](T& self) { return self.invert_options(); });
    c.def(
        "invert_options", [](T& self, const std::string& tpe) { return self.invert_options(tpe); }, "type"_a);
    c.def(
        "apply_inverse",
        [](T& self,
           const V& range,
           V& source,
           const XT::Common::Configuration& opts,
           const XT::Common::Parameter& param) { self.apply_inverse(range, source, opts, param); },
        "range"_a,
        "source"_a,
        "opts"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](T& self, const V& range, V& source, const std::string& tpe, const XT::Common::Parameter& param) {
          self.apply_inverse(range, source, tpe, param);
        },
        "range"_a,
        "source"_a,
        "type"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](T& self, const V& range, V& source, const XT::Common::Parameter& param) {
          self.apply_inverse(range, source, param);
        },
        "range"_a,
        "source"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](T& self, const V& range, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
          return self.apply_inverse(range, opts, param);
        },
        "range"_a,
        "opts"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](T& self, const V& range, const std::string& tpe, const XT::Common::Parameter& param) {
          return self.apply_inverse(range, tpe, param);
        },
        "range"_a,
        "type"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](T& self, const V& range, const XT::Common::Parameter& param) { return self.apply_inverse(range, param); },
        "range"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def("jacobian_options", [](T& self) { return self.jacobian_options(); });
    c.def(
        "jacobian_options", [](T& self, const std::string& tpe) { return self.jacobian_options(tpe); }, "type"_a);
    c.def(
        "jacobian",
        [](T& self,
           const V& source,
           Mop& jacobian_op,
           const XT::Common::Configuration& opts,
           const XT::Common::Parameter& param) { self.jacobian(source, jacobian_op, opts, param); },
        "source"_a,
        "jacobian_op"_a,
        "opts"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "jacobian",
        [](T& self, const V& source, Mop& jacobian_op, const std::string& tpe, const XT::Common::Parameter& param) {
          self.jacobian(source, jacobian_op, tpe, param);
        },
        "source"_a,
        "jacobian_op"_a,
        "type"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "jacobian",
        [](T& self, const V& source, Mop& jacobian_op, const XT::Common::Parameter& param) {
          self.jacobian(source, jacobian_op, param);
        },
        "source"_a,
        "jacobian_op"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "jacobian",
        [](T& self, const V& source, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
          return self.jacobian(source, opts, param);
        },
        "source"_a,
        "opts"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "jacobian",
        [](T& self, const V& source, const std::string& tpe, const XT::Common::Parameter& param) {
          return self.jacobian(source, tpe, param);
        },
        "source"_a,
        "type"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "jacobian",
        [](T& self, const V& source, const XT::Common::Parameter& param) { return self.jacobian(source, param); },
        "source"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    // ... induced norm?

    // operators. These are delicate: we need to mimic the C++ situation, e.g. som operators here, others in
    // ConstLinccomb and Lincomb...
    // * variants from OperatorInterface
    // NOTE: the non-const variants come first on purpose, to be tried first
    c.def(
        "__mul__",
        [](T& self, const F& alpha) { return std::make_unique<decltype(self * alpha)>(self * alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    c.def(
        "__mul__",
        [](const T& self, const F& alpha) { return std::make_unique<decltype(self * alpha)>(self * alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    c.def(
        "__truediv__",
        [](T& self, const F& alpha) { return std::make_unique<decltype(self / alpha)>(self / alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    c.def(
        "__truediv__",
        [](const T& self, const F& alpha) { return std::make_unique<decltype(self / alpha)>(self / alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    c.def(
        "__add__",
        [](T& self, Lop& other) { return std::make_unique<decltype(self + other)>(self + other); },
        "lincomb_op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__add__",
        [](T& self, type& other) { return std::make_unique<decltype(self + other)>(self + other); },
        "op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__add__",
        [](T& self, const V& vec) { return std::make_unique<decltype(self + vec)>(self + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__add__",
        [](const T& self, const CLop& other) { return std::make_unique<decltype(self + other)>(self + other); },
        "const_lincomb_op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__add__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self + other)>(self + other); },
        "op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__add__",
        [](const T& self, const V& vec) { return std::make_unique<decltype(self + vec)>(self + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](T& self, Lop& other) { return std::make_unique<decltype(self - other)>(self - other); },
        "lincomb_op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](T& self, type& other) { return std::make_unique<decltype(self - other)>(self - other); },
        "op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](T& self, const V& vec) { return std::make_unique<decltype(self - vec)>(self - vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](const T& self, const CLop& other) { return std::make_unique<decltype(self - other)>(self - other); },
        "const_lincomb_op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](const T& self, const type& other) { return std::make_unique<decltype(self - other)>(self - other); },
        "op"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__sub__",
        [](const T& self, const V& vec) { return std::make_unique<decltype(self - vec)>(self - vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    // * variants from OperatorInterface with interchanged arguments (at most the combinations from above!)
    c.def(
        "__rmul__",
        [](T& self, const F& alpha) { return std::make_unique<decltype(self * alpha)>(self * alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    c.def(
        "__rmul__",
        [](const T& self, const F& alpha) { return std::make_unique<decltype(self * alpha)>(self * alpha); },
        "scalar"_a,
        py::keep_alive<0, 1>(),
        py::is_operator());
    // __rtruediv__ as in scalar/operators does not make sense
    // __radd__ for other ops does not make sense, uses __add__ of the other op
    c.def(
        "__radd__",
        [](T& self, const V& vec) { return std::make_unique<decltype(self + vec)>(self + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__radd__",
        [](const T& self, const V& vec) { return std::make_unique<decltype(self + vec)>(self + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    // __rsub__ for other ops does not make sense, uses __sub__ of the other op
    c.def(
        "__rsub__",
        [](T& self, const V& vec) { return std::make_unique<decltype(self * -1 + vec)>(self * -1 + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    c.def(
        "__rsub__",
        [](const T& self, const V& vec) { return std::make_unique<decltype(self * -1 + vec)>(self * -1 + vec); },
        "vector"_a,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::is_operator());
    // * additional variants for Python which make sense given OperatorInterface
    c.def(
        "__neg__",
        [](T& self) { return std::make_unique<decltype(self * -1)>(self * -1); },
        py::is_operator(),
        py::keep_alive<0, 1>());
    c.def(
        "__neg__",
        [](const T& self) { return std::make_unique<decltype(self * -1)>(self * -1); },
        py::is_operator(),
        py::keep_alive<0, 1>());
  } // ... addbind_methods(...)

  static std::string class_name(const std::string& matrix_id,
                                const std::string& grid_id,
                                const std::string& layer_id,
                                const std::string& class_id)
  {
    std::string ret = class_id;
    ret += "_" + grid_id;
    if (!layer_id.empty())
      ret += "_" + layer_id;
    ret += "_" + XT::Common::to_string(r_r) + "d_range_space";
    ret += "_" + XT::Common::to_string(s_r) + "d_source_space";
    ret += "_" + matrix_id + "_matrix";
    return ret;
  } // ... class_name(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "operator_interface")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(class_name(matrix_id, grid_id, layer_id, class_id));
    bound_type c(m,
                 ClassName.c_str(),
                 (XT::Common::to_camel_case(class_id) + " (" + grid_id + ", " + matrix_id + " variant)").c_str());

    c.def_property_readonly("linear", &type::linear);
    c.def_property_readonly("source_space",
                            py::cpp_function(&type::source_space, py::return_value_policy::reference_internal));
    c.def_property_readonly("range_space",
                            py::cpp_function(&type::range_space, py::return_value_policy::reference_internal));

    addbind_methods(c);

    return c;
  } // ... bind(...)
}; // class OperatorInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_OPERATORS_INTERFACES_HH
