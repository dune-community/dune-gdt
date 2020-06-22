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

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class G, size_t s = 1, size_t sC = 1, size_t r = s, size_t rC = sC>
class OperatorInterface
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using type = GDT::OperatorInterface<M, GV, s, sC, r, rC>;
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

    c.def("assemble", [](T& self, const bool parallel) { self.assemble(parallel); }, "parallel"_a = false);
    c.def("apply",
          [](T& self, const V& source, V& range, const XT::Common::Parameter& param) {
            self.apply(source, range, param);
          },
          "source"_a,
          "range"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply",
          [](T& self, const V& source, const XT::Common::Parameter& param) { return self.apply(source, param); },
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply2",
          [](T& self, const V& range, const V& source, const XT::Common::Parameter& param) {
            self.apply2(range, source, param);
          },
          "range"_a,
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("invert_options", [](T& self) { return self.invert_options(); });
    c.def("invert_options", [](T& self, const std::string& tpe) { return self.invert_options(tpe); }, "type"_a);
    c.def("apply_inverse",
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
    c.def("apply_inverse",
          [](T& self, const V& range, V& source, const std::string& tpe, const XT::Common::Parameter& param) {
            self.apply_inverse(range, source, tpe, param);
          },
          "range"_a,
          "source"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](T& self, const V& range, V& source, const XT::Common::Parameter& param) {
            self.apply_inverse(range, source, param);
          },
          "range"_a,
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](T& self, const V& range, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
            return self.apply_inverse(range, opts, param);
          },
          "range"_a,
          "opts"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](T& self, const V& range, const std::string& tpe, const XT::Common::Parameter& param) {
            return self.apply_inverse(range, tpe, param);
          },
          "range"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](T& self, const V& range, const XT::Common::Parameter& param) { return self.apply_inverse(range, param); },
          "range"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian_options", [](T& self) { return self.jacobian_options(); });
    c.def("jacobian_options", [](T& self, const std::string& tpe) { return self.jacobian_options(tpe); }, "type"_a);
    c.def("jacobian",
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
    c.def("jacobian",
          [](T& self, const V& source, Mop& jacobian_op, const std::string& tpe, const XT::Common::Parameter& param) {
            self.jacobian(source, jacobian_op, tpe, param);
          },
          "source"_a,
          "jacobian_op"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](T& self, const V& source, Mop& jacobian_op, const XT::Common::Parameter& param) {
            self.jacobian(source, jacobian_op, param);
          },
          "source"_a,
          "jacobian_op"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](T& self, const V& source, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
            return self.jacobian(source, opts, param);
          },
          "source"_a,
          "opts"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](T& self, const V& source, const std::string& tpe, const XT::Common::Parameter& param) {
            return self.jacobian(source, tpe, param);
          },
          "source"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](T& self, const V& source, const XT::Common::Parameter& param) { return self.jacobian(source, param); },
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    // ... induced norm?

    // operators (do we need the const variants? if yes, adjust lincomb bindings accordingly)
    // * from OperatorInterface
    //    virtual ConstLincombOperatorType operator*(const FieldType& alpha) const
    //    c.def("__mul__", [](const type& self, const F& alpha) { return self * alpha; }, "alpha"_a, py::is_operator());
    //    c.def("__rmul__", [](const type& self, const F& alpha) { return self * alpha; }, "alpha"_a,
    //    py::is_operator()); virtual LincombOperatorType operator*(const FieldType& alpha)
    c.def("__mul__", [](type& self, const F& alpha) { return self * alpha; }, "alpha"_a, py::is_operator());
    c.def("__rmul__", [](type& self, const F& alpha) { return self * alpha; }, "alpha"_a, py::is_operator());
    c.def("__imul__", [](type& self, const F& alpha) { return self * alpha; }, "alpha"_a, py::is_operator());
    //    virtual ConstLincombOperatorType operator/(const FieldType& alpha) const
    //    c.def("__truediv__", [](const type& self, const F& alpha) { return self / alpha; }, "alpha"_a,
    //    py::is_operator()); virtual LincombOperatorType operator/(const FieldType& alpha)
    c.def("__truediv__", [](type& self, const F& alpha) { return self / alpha; }, "alpha"_a, py::is_operator());
    c.def("__itruediv__", [](type& self, const F& alpha) { return self / alpha; }, "alpha"_a, py::is_operator());
    //    virtual ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const
    //    c.def("__add__",
    //          [](const type& self, const CLop& other) { return self + other; },
    //          "other_const_lincomb_operator"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual ConstLincombOperatorType operator+(const ThisType& other) const
    //    c.def("__add__",
    //          [](const type& self, const type& other) { return self + other; },
    //          "other_const_operator"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual ConstLincombOperatorType operator+(const VectorType& vector) const
    //    c.def("__add__",
    //          [](const type& self, const V& other) { return self + other; },
    //          "other_const_vector"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator+(LincombOperatorType& other)
    c.def("__add__",
          [](type& self, Lop& other) { return self + other; },
          "other_lincomb_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__iadd__",
          [](type& self, Lop& other) { return self + other; },
          "other_lincomb_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator+(ThisType& other)
    c.def("__add__",
          [](type& self, type& other) { return self + other; },
          "other_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__iadd__",
          [](type& self, type& other) { return self + other; },
          "other_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator+(const VectorType& vector)
    c.def("__add__",
          [](type& self, V& other) { return self + other; },
          "other_vector"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__iadd__",
          [](type& self, V& other) { return self + other; },
          "other_vector"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    //    virtual ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const
    //    c.def("__sub__",
    //          [](const type& self, const CLop& other) { return self - other; },
    //          "other_const_lincomb_operator"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual ConstLincombOperatorType operator-(const ThisType& other) const
    //    c.def("__sub__",
    //          [](const type& self, const type& other) { return self - other; },
    //          "other_const_operator"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual ConstLincombOperatorType operator-(const VectorType& vector) const
    //    c.def("__sub__",
    //          [](const type& self, const V& other) { return self - other; },
    //          "other_const_vector"_a,
    //          py::is_operator(),
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator-(LincombOperatorType& other)
    c.def("__sub__",
          [](type& self, Lop& other) { return self - other; },
          "other_lincomb_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__isub__",
          [](type& self, Lop& other) { return self - other; },
          "other_lincomb_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator-(ThisType& other)
    c.def("__sub__",
          [](type& self, type& other) { return self - other; },
          "other_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__isub__",
          [](type& self, type& other) { return self - other; },
          "other_operator"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    //    virtual LincombOperatorType operator-(const VectorType& vector)
    c.def("__sub__",
          [](type& self, V& other) { return self - other; },
          "other_vector"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("__isub__",
          [](type& self, V& other) { return self - other; },
          "other_vector"_a,
          py::is_operator(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    // others
    //    c.def("__neg__", [](const type& self) { return self * -1; }, py::is_operator(), py::keep_alive<0, 1>());
    c.def("__neg__", [](type& self) { return self * -1; }, py::is_operator(), py::keep_alive<0, 1>());

  } // ... addbind_methods(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "operator_interface",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto matrix_name = XT::LA::bindings::container_name<M>::value();
    //    auto vector_name = XT::LA::bindings::container_name<V>::value();

    std::string class_name = matrix_name + "_" + class_id + "_" + grid_id + "_" + XT::Common::to_string(s);
    if (sC > 1)
      class_name += "x" + XT::Common::to_string(sC);
    class_name += "d_source_space";
    class_name += XT::Common::to_string(r);
    if (rC > 1)
      class_name += "x" + XT::Common::to_string(rC);
    class_name += "d_range_space";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m,
                 ClassName.c_str(),
                 (XT::Common::to_camel_case(class_id) + " (" + grid_id + ", " + matrix_name + " variant)").c_str());

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
