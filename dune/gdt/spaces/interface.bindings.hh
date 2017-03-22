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
//#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>

#include "interface.hh"
#include "cg.hh"
#include "dg.hh"
#include "fv.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <ChooseSpaceBackend backend>
struct backend_name
{
  static_assert(AlwaysFalse<typename internal::backend_dependent_typename<backend>::type>::value,
                "Please add a specialization for this backend!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct backend_name<ChooseSpaceBackend::fem>
{
  static std::string value()
  {
    return "fem";
  }
};

template <>
struct backend_name<ChooseSpaceBackend::gdt>
{
  static std::string value()
  {
    return "gdt";
  }
};

template <>
struct backend_name<ChooseSpaceBackend::pdelab>
{
  static std::string value()
  {
    return "pdelab";
  }
};


namespace internal {


template <class G, XT::Grid::Layers layer, ChooseSpaceBackend backend, size_t r, size_t rC>
struct space_name_base
{
  static std::string value()
  {
    using XT::Common::to_string;
    return XT::Grid::bindings::grid_name<G>::value() + "_" + XT::Grid::bindings::layer_name<layer>::value() + "_to_"
           + to_string(r) + "x" + to_string(rC) + "_" + backend_name<backend>::value();
  }
};


} // namespace internal


template <class P>
struct space_name
{
  static_assert(AlwaysFalse<P>::value, "Please add a specialization for this space provider!");

  static std::string value()
  {
    return "";
  }
};

template <class G, XT::Grid::Layers layer, ChooseSpaceBackend backend, int p, size_t r, size_t rC>
struct space_name<CgSpaceProvider<G, layer, backend, p, double, r, rC>>
{
  static std::string value()
  {
    return std::string("cg_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};

template <class G, XT::Grid::Layers layer, ChooseSpaceBackend backend, int p, size_t r, size_t rC>
struct space_name<DgSpaceProvider<G, layer, backend, p, double, r, rC>>
{
  static std::string value()
  {
    return std::string("dg_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};

template <class G, XT::Grid::Layers layer, ChooseSpaceBackend backend, size_t r, size_t rC>
struct space_name<FvSpaceProvider<G, layer, backend, double, r, rC>>
{
  static std::string value()
  {
    return std::string("fv_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_space";
  }
};


template <class SP>
class SpaceInterface
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");

public:
  typedef S type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = space_name<SP>::value();

    bound_type c(m, ClassName.c_str(), ClassName.c_str(), py::metaclass());

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

//#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
