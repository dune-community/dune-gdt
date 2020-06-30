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

public:
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "const_lincomb_operator",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value())
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
    c.def("add",
          [](type& self, const typename type::OperatorType& op, const typename type::FieldType& coeff) {
            self.add(op, coeff);
          },
          "op"_a,
          "coeff"_a = 1.);
    c.def("add",
          [](type& self, const type& op, const typename type::FieldType& coeff) { self.add(op, coeff); },
          "const_lincomb_op"_a,
          "coeff"_a = 1.);
    c.def_property_readonly("num_ops", &type::num_ops);
    c.def("op", [](type& self, const size_t ii) { return self.op(ii); }, "index"_a, py::keep_alive<1, 0>());
    c.def("coeff", [](type& self, const size_t ii) { return self.coeff(ii); }, "index"_a);
    // operators
    // ...

    // factory
    //    m.def(XT::Common::to_camel_case(class_id).c_str(),
    //          [](const SS& source_space, const RS& range_space, const MT&) {
    //            return type(source_space, range_space);
    //          },
    //          "source_space"_a,
    //          "range_space"_a,
    //          "matrix_type"_a,
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    m.def(XT::Common::to_camel_case(class_id).c_str(),
          [](const SS& source_space, const RS& range_space, const M&) { return type(source_space, range_space); },
          "source_space"_a,
          "range_space"_a,
          "unsused_matrix_to_select_type"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    m.def(XT::Common::to_camel_case(matrix_name + "_" + class_id).c_str(),
          [](const SS& source_space, const RS& range_space) { return type(source_space, range_space); },
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

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
  using Mop = typename type::MatrixOperatorType;

public:
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "lincomb_operator",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value())
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
    c.def("add",
          [](type& self, const typename type::OperatorType& op, const typename type::FieldType& coeff) {
            self.add(op, coeff);
          },
          "op"_a,
          "coeff"_a = 1.);
    c.def(
        "add",
        [](type& self, typename type::OperatorType& op, const typename type::FieldType& coeff) { self.add(op, coeff); },
        "op"_a,
        "coeff"_a = 1.);
    c.def("add",
          [](type& self, const type& op, const typename type::FieldType& coeff) { self.add(op, coeff); },
          "const_lincomb_op"_a,
          "coeff"_a = 1.);
    c.def("add",
          [](type& self, type& op, const typename type::FieldType& coeff) { self.add(op, coeff); },
          "const_lincomb_op"_a,
          "coeff"_a = 1.);
    c.def_property_readonly("num_ops", &type::num_ops);
    c.def("op", [](type& self, const size_t ii) { return self.op(ii); }, "index"_a, py::keep_alive<1, 0>());
    c.def("coeff", [](type& self, const size_t ii) { return self.coeff(ii); }, "index"_a);
    // operators
    // ...

    // factory
    //    m.def(XT::Common::to_camel_case(class_id).c_str(),
    //          [](const SS& source_space, const RS& range_space, const MT&) {
    //            return type(source_space, range_space);
    //          },
    //          "source_space"_a,
    //          "range_space"_a,
    //          "matrix_type"_a,
    //          py::keep_alive<0, 1>(),
    //          py::keep_alive<0, 2>());
    m.def(XT::Common::to_camel_case(class_id).c_str(),
          [](const SS& source_space, const RS& range_space, const M&) { return type(source_space, range_space); },
          "source_space"_a,
          "range_space"_a,
          "unsused_matrix_to_select_type"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    m.def(XT::Common::to_camel_case(matrix_name + "_" + class_id).c_str(),
          [](const SS& source_space, const RS& range_space) { return type(source_space, range_space); },
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class LincombOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_OPERATORS_LINCOMB_HH
