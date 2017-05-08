// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_BINDINGS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include "assembler.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class LocalCouplingTwoFormType>
class LocalCouplingTwoFormAssembler
{
  static_assert(is_local_coupling_twoform<LocalCouplingTwoFormType>::value, "");

public:
  typedef GDT::LocalCouplingTwoFormAssembler<LocalCouplingTwoFormType> type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& local_coupling_two_form_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("local_coupling_" + local_coupling_two_form_name + "_assembler");

    bound_type c(m, ClassName.c_str());

    m.def("make_local_coupling_two_form_assembler",
          [](const LocalCouplingTwoFormType& local_coupling_two_form) { return type(local_coupling_two_form); },
          "local_coupling_two_form"_a);

    return c;
  } // ... bind(...)
}; // class LocalCouplingTwoFormAssembler


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_LOCAL_ASSEMBLER_BINDINGS_HH
