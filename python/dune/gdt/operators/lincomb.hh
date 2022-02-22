// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_OPERATORS_LINCOMB_HH
#define PYTHON_DUNE_GDT_OPERATORS_LINCOMB_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/operators/lincomb.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class GV, size_t s = 1, size_t r = s>
class ConstLincombOperator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  using type = GDT::ConstLincombOperator<GV, s, 1, r, 1>; // M
  using base_type = GDT::OperatorInterface<GV, s, 1, r, 1>; // M
  using V = typename type::VectorType;
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using Mop = typename type::MatrixOperatorType;
  using F = typename type::FieldType;

public:
  using bound_type = pybind11::class_<type, base_type>;

  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    // methods from operator base, to allow for overloads
    bindings::OperatorInterface<M, GV, s, r>::addbind_methods(c); // M

    // our methods
    c.def(
        "add",
        [](T& self, const typename type::OperatorType& op, const F& coeff) { self.add(op, coeff); },
        "op"_a,
        "coeff"_a = 1.,
        py::keep_alive<1, 2>());
    c.def(
        "add",
        [](T& self, const type& op, const F& coeff) { self.add(op, coeff); },
        "const_lincomb_op"_a,
        "coeff"_a = 1.,
        py::keep_alive<1, 2>());
    c.def(
        "op", (const typename T::OperatorType& (T::*)(const size_t) const) & T::op, "index"_a, py::keep_alive<0, 1>());
    c.def(
        "coeff", [](T& self, const size_t ii) { return self.coeff(ii); }, "index"_a);

    // our operators
    // (function ptr signature required for the right return type)
    c.def("__imul__", (type & (type::*)(const F&)) & type::operator*=, py::is_operator());
    c.def("__itruediv__", (type & (type::*)(const F&)) & type::operator/=, py::is_operator());
    c.def(
        "__iadd__", (type & (type::*)(const base_type&)) & type::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__iadd__", (type & (type::*)(const type&)) & type::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def(
        "__isub__", (type & (type::*)(const base_type&)) & type::operator-=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__isub__", (type & (type::*)(const type&)) & type::operator-=, py::is_operator(), py::keep_alive<1, 2>());
  } // ... addbind_methods(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "const_lincomb_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        bindings::OperatorInterface<M, GV, s, r>::class_name(matrix_id, grid_id, layer_id, class_id));
    bound_type c(m,
                 ClassName.c_str(),
                 (XT::Common::to_camel_case(class_id) + " (" + grid_id + ", " + matrix_id + " variant)").c_str());
    c.def(py::init<const SS&, const RS&>(),
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());
    c.def_property_readonly("num_ops", &type::num_ops);

    addbind_methods(c);

    return c;
  } // ... bind(...)
}; // class ConstLincombOperator


template <class M, class GV, size_t s = 1, size_t r = s>
class LincombOperator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  using type = GDT::LincombOperator<GV, s, 1, r, 1>; // M
  using base_type = GDT::ConstLincombOperator<GV, s, 1, r, 1>; // M
  using V = typename type::VectorType;
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using F = typename type::FieldType;
  using Mop = typename type::MatrixOperatorType;

public:
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "lincomb_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        bindings::OperatorInterface<M, GV, s, r>::class_name(matrix_id, grid_id, layer_id, class_id)); // M
    bound_type c(m,
                 ClassName.c_str(),
                 (XT::Common::to_camel_case(class_id) + " (" + grid_id + ", " + matrix_id + " variant)").c_str());
    c.def(py::init<const SS&, const RS&>(),
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    // methods from base, to allow for overloads
    bindings::ConstLincombOperator<M, GV, s, r>::addbind_methods(c);

    // our methods
    c.def(
        "add",
        [](type& self, typename type::OperatorType& op, const F& coeff) { self.add(op, coeff); },
        "op"_a,
        "coeff"_a = 1.,
        py::keep_alive<1, 2>());
    c.def(
        "add",
        [](type& self, type& op, const F& coeff) { self.add(op, coeff); },
        "const_lincomb_op"_a,
        "coeff"_a = 1.,
        py::keep_alive<1, 2>());
    c.def("op", (typename type::OperatorType & (type::*)(const size_t)) & type::op, "index"_a, py::keep_alive<0, 1>());
    c.def("assemble",
          (typename type::OperatorType & (type::*)(const bool)) & type::assemble,
          "parallel"_a = false,
          py::keep_alive<0, 1>());

    // our operators, only those that are not yet present in OperatorInterface or ConstLincombOperator
    // (function ptr signature required for the right return type)
    c.def("__iadd__", (type & (type::*)(type&)) & type::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__isub__", (type & (type::*)(type&)) & type::operator-=, py::is_operator(), py::keep_alive<1, 2>());

    return c;
  } // ... bind(...)
}; // class LincombOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_OPERATORS_LINCOMB_HH
