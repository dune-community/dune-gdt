// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/grids.bindings.hh>

#include <dune/gdt/spaces.bindings.hh>

#include "system.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class T /*, class GV = typename T::GridViewType, class A = T*/>
class SystemAssembler
{
  static_assert(is_space<T>::value, "");

public:
  typedef GDT::SystemAssembler<T /*, G, A*/> type;
  typedef pybind11::class_<type> bound_type;

private:
  typedef typename type::TestSpaceType TestSpaceType;
  typedef typename type::GridViewType GridViewType;
  typedef typename type::AnsatzSpaceType AnsatzSpaceType;

  template <bool do_bind = (std::is_same<TestSpaceType, AnsatzSpaceType>::value
                            && std::is_same<GridViewType, typename TestSpaceType::GridViewType>::value),
            bool anything = true>
  struct addbind_ctor_single
  {
    void operator()(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      c.def(py::init<T>(),
            "space"_a,
            "Uses given space as test and ansatz space, and the grid view of the given space as grid view.",
            py::keep_alive<1, 2>());
    }
  }; // struct addbind_ctor_single

  template <bool anything>
  struct addbind_ctor_single<false, anything>
  {
    void operator()(bound_type& /*c*/)
    {
    }
  }; // struct addbind_ctor_single

public:
  static bound_type bind(pybind11::module& m, const std::string& space_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(
        m, std::string("SystemAssembler__" + space_id).c_str(), std::string("SystemAssembler__" + space_id).c_str());
    addbind_ctor_single<>()(c);

    c.def("append", [](type& self, type& other) { self.append(other); }, "system_assembler"_a, py::keep_alive<1, 2>());
    c.def("append",
          [](type& self, XT::Grid::Walker<GridViewType>& other) { self.append(other); },
          "grid_walker"_a,
          py::keep_alive<1, 2>());
    c.def("assemble", [](type& self, const bool use_tbb) { self.assemble(use_tbb); }, "use_tbb"_a = false);

    m.def("make_system_assembler",
          [](const TestSpaceType& space) { return new type(space); },
          "space"_a,
          py::keep_alive<0, 1>());
    return c;
  } // ... bind(...)
}; // class SystemAssembler


#define DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(_prefix, _GRID)                                                             \
  _prefix class SystemAssembler<FV_SPACE(_GRID, leaf, gdt, 1, 1)>;                                                     \
  _prefix class SystemAssembler<FV_SPACE(_GRID, level, gdt, 1, 1)>

#define DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(_prefix, _GRID)                                                             \
  _prefix class SystemAssembler<CG_SPACE(_GRID, leaf, fem, 1, 1, 1)>;                                                  \
  _prefix class SystemAssembler<CG_SPACE(_GRID, level, fem, 1, 1, 1)>;                                                 \
  _prefix class SystemAssembler<DG_SPACE(_GRID, leaf, fem, 1, 1, 1)>;                                                  \
  _prefix class SystemAssembler<DG_SPACE(_GRID, level, fem, 1, 1, 1)>;                                                 \
  _prefix class SystemAssembler<DG_SPACE(_GRID, leaf, fem, 2, 1, 1)>;                                                  \
  _prefix class SystemAssembler<DG_SPACE(_GRID, level, fem, 2, 1, 1)>

#define DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(_prefix, _GRID)                                                          \
  _prefix class SystemAssembler<CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1)>;                                               \
  _prefix class SystemAssembler<CG_SPACE(_GRID, level, pdelab, 1, 1, 1)>


// these lines have to match the corresponding ones in the .cc source file
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#if HAVE_DUNE_FEM
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_DUNE_ALUGRID
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(extern template, ALU_2D_SIMPLEX_CONFORMING);
#if HAVE_DUNE_FEM
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif // HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
