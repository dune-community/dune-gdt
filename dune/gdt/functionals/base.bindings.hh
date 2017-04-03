// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_FUNCTIONALS_BASE_BINDINGS_H
#define DUNE_GDT_FUNCTIONALS_BASE_BINDINGS_H
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class FunctionalType>
class VectorFunctionalBase
{
public:
  typedef FunctionalType type;
  typedef GDT::SystemAssembler<typename FunctionalType::SpaceType, typename FunctionalType::GridLayerType> BaseType;
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& class_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    typedef typename type::SpaceType S;
    typedef typename type::VectorType V;

    bound_type c(m, std::string(class_id).c_str(), std::string(class_id).c_str());

    c.def("vector", [](type& self) { return self.vector(); });
    c.def("space", [](type& self) { return self.space(); });
    c.def("apply", [](type& self, const V& source) { return self.apply(source); }, "source"_a);
    c.def(
        "apply", [](type& self, const ConstDiscreteFunction<S, V>& source) { return self.apply(source); }, "source"_a);

    return c;
  } // ... bind(...)
}; // class VectorFunctionalBase


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_FUNCTIONALS_BASE_BINDINGS_H
