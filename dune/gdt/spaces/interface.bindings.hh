// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
#define DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class S>
class SpaceInterface
{
  static_assert(is_space<S>::value, "");

public:
  typedef S type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(m, id.c_str(), id.c_str(), py::metaclass());

    c.def_property_readonly_static("dimDomain", [](const type& /*self*/) { return S::dimDomain; });
    c.def_property_readonly_static("dimRange", [](const type& /*self*/) { return S::dimRange; });
    c.def_property_readonly_static("dimRangeCols", [](const type& /*self*/) { return S::dimRangeCols; });
    c.def_property_readonly_static("polOrder", [](const type& /*self*/) { return S::polOrder; });

    c.def("size", [](const type& self) { return self.mapper().size(); });
    c.def("visualize",
          [](const type& self, const std::string& filename) { self.visualize(filename); },
          "filename"_a = "");
    c.def("compute_pattern", [](const type& self) { return self.compute_pattern(); });
    c.def("compute_volume_pattern", [](const type& self) { return self.compute_volume_pattern(); });
    c.def("compute_face_pattern", [](const type& self) { return self.compute_face_pattern(); });
    c.def("compute_face_and_volume_pattern", [](const type& self) { return self.compute_face_and_volume_pattern(); });

    return c;
  } // ... bind(...)
}; // class SpaceInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
