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
  using Mop = typename type::MatrixOperatorType;

public:
  using bound_type = pybind11::class_<type>;

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
    c.def_property_readonly("source_space", &type::source_space);
    c.def_property_readonly("range_space", &type::range_space);
    c.def("assemble", [](type& self, const bool parallel) { self.assemble(parallel); }, "parallel"_a = false);
    c.def("apply",
          [](type& self, const V& source, V& range, const XT::Common::Parameter& param) {
            self.apply(source, range, param);
          },
          "source"_a,
          "range"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply",
          [](type& self, const V& source, const XT::Common::Parameter& param) { return self.apply(source, param); },
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply2",
          [](type& self, const V& range, const V& source, const XT::Common::Parameter& param) {
            self.apply2(range, source, param);
          },
          "range"_a,
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("invert_options", [](type& self) { return self.invert_options(); });
    c.def("invert_options", [](type& self, const std::string& tpe) { return self.invert_options(tpe); }, "type"_a);
    c.def("apply_inverse",
          [](type& self,
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
          [](type& self, const V& range, V& source, const std::string& tpe, const XT::Common::Parameter& param) {
            self.apply_inverse(range, source, tpe, param);
          },
          "range"_a,
          "source"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](type& self, const V& range, V& source, const XT::Common::Parameter& param) {
            self.apply_inverse(range, source, param);
          },
          "range"_a,
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](type& self, const V& range, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
            return self.apply_inverse(range, opts, param);
          },
          "range"_a,
          "opts"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("apply_inverse",
          [](type& self, const V& range, const std::string& tpe, const XT::Common::Parameter& param) {
            return self.apply_inverse(range, tpe, param);
          },
          "range"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def(
        "apply_inverse",
        [](type& self, const V& range, const XT::Common::Parameter& param) { return self.apply_inverse(range, param); },
        "range"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def("jacobian_options", [](type& self) { return self.jacobian_options(); });
    c.def("jacobian_options", [](type& self, const std::string& tpe) { return self.jacobian_options(tpe); }, "type"_a);
    c.def("jacobian",
          [](type& self,
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
        [](type& self, const V& source, Mop& jacobian_op, const std::string& tpe, const XT::Common::Parameter& param) {
          self.jacobian(source, jacobian_op, tpe, param);
        },
        "source"_a,
        "jacobian_op"_a,
        "type"_a,
        "param"_a = XT::Common::Parameter(),
        py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](type& self, const V& source, Mop& jacobian_op, const XT::Common::Parameter& param) {
            self.jacobian(source, jacobian_op, param);
          },
          "source"_a,
          "jacobian_op"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](type& self, const V& source, const XT::Common::Configuration& opts, const XT::Common::Parameter& param) {
            return self.jacobian(source, opts, param);
          },
          "source"_a,
          "opts"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](type& self, const V& source, const std::string& tpe, const XT::Common::Parameter& param) {
            return self.jacobian(source, tpe, param);
          },
          "source"_a,
          "type"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    c.def("jacobian",
          [](type& self, const V& source, const XT::Common::Parameter& param) { return self.jacobian(source, param); },
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
    // ... induced norm?

    return c;
  } // ... bind(...)
}; // class OperatorInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class M, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct OperatorInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::OperatorInterface<M, G>::bind(m);
    if (d > 1) {
      Dune::GDT::bindings::OperatorInterface<M, G, d, 1, 1, 1>::bind(m);
      Dune::GDT::bindings::OperatorInterface<M, G, 1, 1, d, 1>::bind(m);
      Dune::GDT::bindings::OperatorInterface<M, G, d, 1, d, 1>::bind(m);
    }
    // add your extra dimensions here
    // ...
    OperatorInterface_for_all_grids<M, typename GridTypes::tail_type>::bind(m);
  }
};

template <class M>
struct OperatorInterface_for_all_grids<M, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


#endif // PYTHON_DUNE_GDT_OPERATORS_INTERFACES_HH
