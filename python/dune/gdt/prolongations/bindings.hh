// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PROLONGATIONS_BINDINGS_HH
#define DUNE_GDT_PROLONGATIONS_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <python/dune/gdt/discretefunction/bindings.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/prolongations.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class SS, class SV, class RS, class RV>
class prolong
{
  static_assert(is_space<SS>::value, "");
  static_assert(is_space<RS>::value, "");
  static_assert(XT::LA::is_vector<SV>::value, "");
  static_assert(XT::LA::is_vector<RV>::value, "");

public:
  static void bind(pybind11::module& m)
  {
    using namespace pybind11::literals;

    m.def("prolong",
          [](const GDT::ConstDiscreteFunction<SS, SV>& source,
             GDT::DiscreteFunction<RS, RV>& range,
             const size_t over_integrate) { Dune::GDT::prolong(source, range, over_integrate); },
          "source"_a,
          "range"_a,
          "over_integrate"_a = 0);
  } // ... bind(...)
}; // class prolong


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_PROLONGATIONS_BINDINGS_HH
