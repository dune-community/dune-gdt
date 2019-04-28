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

#ifndef PYTHON_DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
#define PYTHON_DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH

//#if  HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <python/dune/xt/la/container.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/functions/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class GV, size_t r = 1, size_t rC = 1, class R = double>
class DiscreteFunction
{
public:
  using type = GDT::DiscreteFunction<V, GV, r, rC, R>;
  using base_type = XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using G = typename GV::Grid;
    using S = typename type::SpaceType;

    // base_type
    const auto grid_name = XT::Grid::bindings::grid_name<G>::value();
    if (std::is_same<R, double>::value)
      XT::Common::bindings::try_register(
          m, [&](auto& module) { XT::Functions::bind_GridFunctionInterface<G, r, rC>(module, grid_name); });

    // type
    std::string ClassName = "DiscreteFunction_" + grid_name;
    if (r > 1)
      ClassName += "_to" + XT::Common::to_string(r);
    if (rC > 1)
      ClassName += "x" + XT::Common::to_string(rC);
    ClassName += "_" + XT::LA::bindings::container_name<V>::value();
    ClassName = XT::Common::to_camel_case(ClassName);
    const std::string default_name = "dune.gdt.discretefunction";
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const S&, V&, const std::string&>(),
          "space"_a,
          "vector"_a,
          "name"_a = default_name,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def(py::init<const S&, const std::string&>(), "space"_a, "name"_a = default_name, py::keep_alive<1, 2>());
    c.def_property_readonly("space", &type::space);
    c.def_property("dof_vector",
                   [](const type& self) { return self.dofs().vector(); },
                   [](type& self, const V& vec) {
                     DUNE_THROW_IF(vec.size() != self.dofs().vector().size(),
                                   Exceptions::discrete_function_error,
                                   "vec.size() = " << vec.size() << "\n   self.dofs().vector().size() = "
                                                   << self.dofs().vector().size());
                     self.dofs().vector() = vec;
                   },
                   py::return_value_policy::reference_internal);
    c.def_property_readonly("name", &type::name);
    c.def("visualize",
          [](type& self, const std::string& filename) { return self.visualize(filename, VTK::appendedraw); },
          "filename"_a);

    // factory
    const std::string factory_name = "make_discrete_function";
    m.def(
        factory_name.c_str(),
        [](const S& space, V& vector, const std::string& name) { return make_discrete_function(space, vector, name); },
        "space"_a,
        "vector"_a,
        "name"_a = default_name,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
    m.def(factory_name.c_str(),
          [](const S& space, const std::string& name) { return make_discrete_function<V>(space, name); },
          "space"_a,
          "name"_a = default_name,
          py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)
}; // class DiscreteFunction


} // namespace bindings
} // namespace GDT
} // namespace Dune


//#endif // HAVE_DUNE_PYBINDXI

#endif // PYTHON_DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
