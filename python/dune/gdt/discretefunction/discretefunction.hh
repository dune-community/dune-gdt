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

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/functions/interfaces/grid-function.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class MatchingVectorTag, class GV, size_t r = 1, size_t rC = 1, class R = double>
class DiscreteFunction
{
  using type = GDT::DiscreteFunction<V, GV, r, rC, R>;
  using base_type = XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>;
  using G = typename GV::Grid;
  using S = typename type::SpaceType;

public:
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "discrete_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    if (r > 1)
      class_name += "_to_" + XT::Common::to_string(r) + "d";
    if (rC > 1)
      class_name += "x" + XT::Common::to_string(rC) + "d";
    class_name += "_" + XT::LA::bindings::container_name<V>::value();
    const auto ClassName = XT::Common::to_camel_case(class_name);
    const std::string default_name = "dune.gdt.discretefunction";
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def(py::init<const S&, V&, const std::string&>(),
          "space"_a,
          "vector"_a,
          "name"_a = default_name,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def(py::init<const S&, const std::string&>(), "space"_a, "name"_a = default_name, py::keep_alive<1, 2>());

    c.def_property_readonly("dim_domain", [](const type&) { return size_t(base_type::d); });
    if (rC == 1)
      c.def_property_readonly("dim_range", [](const type&) { return size_t(r); });
    else
      c.def_property_readonly("dim_range", [](const type&) { return std::make_pair(size_t(r), size_t(rC)); });
    c.def_property_readonly("space", &type::space);
    c.def_property("dofs", // doing this so complicated to get an actual reference instead of a copy
                   (typename type::DofVectorType & (type::*)()) & type::dofs,
                   (typename type::DofVectorType & (type::*)()) & type::dofs);
    c.def_property_readonly("name", &type::name);

    c.def(
        "visualize",
        [](type& self, const std::string& filename) { return self.visualize(filename, VTK::appendedraw); },
        "filename"_a);

    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const S& space, V& vector, const std::string& name) { return make_discrete_function(space, vector, name); },
        "space"_a,
        "vector"_a,
        "name"_a = default_name,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const S& space, const MatchingVectorTag&, const std::string& name) {
          return make_discrete_function(space, name);
        },
        "space"_a,
        "vector_type"_a,
        "name"_a = default_name,
        py::keep_alive<0, 1>());
    if (std::is_same<MatchingVectorTag, XT::LA::bindings::Istl>::value)
      m.def(
          XT::Common::to_camel_case(class_id).c_str(),
          [](const S& space, const std::string& name, const MatchingVectorTag&) {
            return make_discrete_function(space, name);
          },
          "space"_a,
          "name"_a = default_name,
          "vector_type"_a = XT::LA::bindings::Istl(),
          py::keep_alive<0, 1>());
    return c;
  } // ... bind(...)
}; // class DiscreteFunction


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
