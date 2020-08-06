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
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/operators/lincomb.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class G, size_t s = 1, size_t sC = 1, size_t r = s, size_t rC = sC>
class ConstLincombOperator
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using type = GDT::ConstLincombOperator<M, GV, s, sC, r, rC>;
  using base_type = GDT::OperatorInterface<M, GV, s, sC, r, rC>;
  using V = typename type::VectorType;
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using Mop = typename type::MatrixOperatorType;
  using F = F;

public:
  using bound_type = pybind11::class_<type, base_type>;

  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    // methods from operator base, to allow for overloads
    bindings::OperatorInterface<M, G, s, sC, r, rC>::addbind_methods(c);

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

    // our operators, only those that are not yet present in OperatorInterface
    // (function ptr signature required for the right return type)
    c.def("__imul__", (T & (T::*)(const F&)) & T::operator*=, py::is_operator());
    c.def("__itruediv__", (T & (T::*)(const F&)) & T::operator/=, py::is_operator());
    c.def("__iadd__", (T & (T::*)(const base_type&)) & T::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__iadd__", (T & (T::*)(const type&)) & T::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__isub__", (T & (T::*)(const base_type&)) & T::operator-=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__isub__", (T & (T::*)(const type&)) & T::operator-=, py::is_operator(), py::keep_alive<1, 2>());
  } // ... addbind_methods(...)

  static bound_type
  bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "const_lincomb_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto matrix_name = XT::LA::bindings::container_name<M>::value();

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


template <class M, class G, size_t s = 1, size_t sC = 1, size_t r = s, size_t rC = sC>
class LincombOperator
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using type = GDT::LincombOperator<M, GV, s, sC, r, rC>;
  using base_type = GDT::ConstLincombOperator<M, GV, s, sC, r, rC>;
  using V = typename type::VectorType;
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using F = F;
  using Mop = typename type::MatrixOperatorType;

public:
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type
  bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "lincomb_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto matrix_name = XT::LA::bindings::container_name<M>::value();

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
    c.def(py::init<const SS&, const RS&>(),
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    // methods from base, to allow for overloads
    bindings::ConstLincombOperator<M, G, s, sC, r, rC>::addbind_methods(c);

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
          "parallel"_a = False,
          py::keep_alive<0, 1>());

    // our operators, only those that are not yet present in OperatorInterface or ConstLincombOperator
    // (function ptr signature required for the right return type)
    c.def("__iadd__", (type & (type::*)(const type&)) & type::operator+=, py::is_operator(), py::keep_alive<1, 2>());
    c.def("__isub__", (type & (type::*)(const type&)) & type::operator-=, py::is_operator(), py::keep_alive<1, 2>());

    return c;
  } // ... bind(...)
}; // class LincombOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_OPERATORS_LINCOMB_HH
