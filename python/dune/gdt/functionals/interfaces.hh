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

#ifndef PYTHON_DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define PYTHON_DUNE_GDT_FUNCTIONALS_INTERFACES_HH

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

#include <dune/gdt/functionals/interfaces.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class G, size_t s = 1, size_t sC = 1, class F = double>
class FunctionalInterface
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using type = GDT::FunctionalInterface<V, GV, s, sC, F>;
  using SS = typename type::SourceSpaceType;

public:
  using bound_type = pybind11::class_<type>;

  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def("assemble", [](type& self, const bool parallel) { self.assemble(parallel); }, "parallel"_a = false);
    c.def("apply",
          [](type& self, const V& source, const XT::Common::Parameter& param) { return self.apply(source, param); },
          "source"_a,
          "param"_a = XT::Common::Parameter(),
          py::call_guard<py::gil_scoped_release>());
  } // ... addbind_methods(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "functional_interface")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto vector_name = XT::LA::bindings::container_name<V>::value();

    std::string class_name = vector_name + "_" + class_id + "_" + grid_id + "_" + XT::Common::to_string(s);
    if (sC > 1)
      class_name += "x" + XT::Common::to_string(sC);
    class_name += "d_source_space";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m,
                 ClassName.c_str(),
                 (XT::Common::to_camel_case(class_id) + " (" + grid_id + ", " + vector_name + " variant)").c_str());

    c.def_property_readonly("linear", &type::linear);
    c.def_property_readonly("source_space", &type::source_space);

    return c;
  } // ... bind(...)
}; // class FunctionalInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct FunctionalInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::FunctionalInterface;

    FunctionalInterface<V, G>::bind(m);
    if (d > 1)
      FunctionalInterface<V, G, d>::bind(m);
    // add your extra dimensions here
    // ...
    FunctionalInterface_for_all_grids<V, typename GridTypes::tail_type>::bind(m);
  }
};

template <class V>
struct FunctionalInterface_for_all_grids<V, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


#endif // PYTHON_DUNE_GDT_FUNCTIONALS_INTERFACES_HH
