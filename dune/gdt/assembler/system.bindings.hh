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
#include <dune/gdt/spaces/constraints.bindings.hh>

#include "system.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class TP /*, class GL = typename TP::type::GridLayerType, class AP = TP*/>
class SystemAssembler
{
  typedef typename TP::type T;
  static_assert(is_space<T>::value, "");

public:
  typedef GDT::SystemAssembler<T /*, GL, A*/> type;
  typedef pybind11::class_<type> bound_type;

private:
  typedef typename type::TestSpaceType TestSpaceType;
  typedef typename type::GridLayerType GridLayerType;
  typedef typename type::AnsatzSpaceType AnsatzSpaceType;

  template <bool do_bind = (std::is_same<TestSpaceType, AnsatzSpaceType>::value
                            && std::is_same<GridLayerType, typename TestSpaceType::GridLayerType>::value),
            bool anything = true>
  struct addbind_ctor_single
  {
    void operator()(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      c.def(py::init<T>(),
            "space"_a,
            "Uses given space as test and ansatz space, and the grid layer of the given space as grid layer.",
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
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(std::string("system_assembler_") + space_name<TP>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    addbind_ctor_single<>()(c);

    c.def("append", [](type& self, type& other) { self.append(other); }, "system_assembler"_a, py::keep_alive<1, 2>());
    c.def("append",
          [](type& self, XT::Grid::Walker<GridLayerType>& other) { self.append(other); },
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


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#if HAVE_DUNE_FEM

// * fem backend
//   - cg spaces

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG(_m, _GRID, _layer, _p, _r, _rC)                                         \
  auto system_assembler_##_GRID##_##_layer##_fem_cg_##_p##_double_##_r##_##_rC =                                       \
      Dune::GDT::bindings::SystemAssembler<Dune::GDT::CgSpaceProvider<_GRID,                                           \
                                                                      Dune::XT::Grid::Layers::_layer,                  \
                                                                      Dune::GDT::Backends::fem,                        \
                                                                      _p,                                              \
                                                                      double,                                          \
                                                                      _r,                                              \
                                                                      _rC>>::bind(_m);                                 \
  Dune::GDT::bindings::                                                                                                \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::part,                          \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::addbind(system_assembler_##_GRID##_##_layer##_fem_cg_##_p##_double_##_r##_##_rC)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALBERTA(_m, _layer, _p)                                               \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG(_m, ALBERTA_2D, _layer, _p, 1, 1)
//#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALBERTA(_m, _p, _layer)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALU(_m, _layer, _p)                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALU(_m, _layer, _p)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_UG(_m, _layer, _p)                                                    \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG(_m, UG_2D, _layer, _p, 1, 1)
//#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_UG(_m, _layer, _p)
//#endif

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_YASP(_m, _layer, _p)                                                    \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1)

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALL(_m, _layer, _p)                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALBERTA(_m, _layer, _p);                                                      \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALU(_m, _layer, _p);                                                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_UG(_m, _layer, _p);                                                           \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_YASP(_m, _layer, _p)

//   - dg spaces

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG(_m, _GRID, _layer, _p, _r, _rC)                                         \
  auto system_assembler_##_GRID##_##_layer##_fem_dg_##_p##_double_##_r##_##_rC =                                       \
      Dune::GDT::bindings::SystemAssembler<Dune::GDT::DgSpaceProvider<_GRID,                                           \
                                                                      Dune::XT::Grid::Layers::_layer,                  \
                                                                      Dune::GDT::Backends::fem,                        \
                                                                      _p,                                              \
                                                                      double,                                          \
                                                                      _r,                                              \
                                                                      _rC>>::bind(_m);                                 \
  Dune::GDT::bindings::                                                                                                \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::part,                          \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::addbind(system_assembler_##_GRID##_##_layer##_fem_dg_##_p##_double_##_r##_##_rC)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALBERTA(_m, _layer, _p)                                               \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG(_m, ALBERTA_2D, _layer, _p, 1, 1)
//#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALBERTA(_m, _p, _layer)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALU(_m, _layer, _p)                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALU(_m, _layer, _p)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG // <- does not work
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_UG(_m, _layer, _p)
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG(_m, UG_2D, _layer, _p, 1, 1)
//#else
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_UG(_m, _layer, _p)
//#endif

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_YASP(_m, _layer, _p)                                                    \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1)

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, _layer, _p)                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALBERTA(_m, _layer, _p);                                                      \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALU(_m, _layer, _p);                                                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_YASP(_m, _layer, _p)
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_UG(_m, _layer, _p); // <- does not work

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(_m)                                                                        \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALL(_m, dd_subdomain, 1);                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALL(_m, leaf, 1);                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_CG_ALL(_m, level, 1);                                                            \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, dd_subdomain, 1);                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, leaf, 1);                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, level, 1);                                                            \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, dd_subdomain, 2);                                                     \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, leaf, 2);                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM_DG_ALL(_m, level, 2)

#else // HAVE_DUNE_FEM
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(_m)
#endif

//#if HAVE_DUNE_PDELAB

//// * pdelab backend
////   - cg spaces

//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG(_m, _GRID, _layer, _p, _r, _rC)                                    \
//  auto system_assembler_##_GRID##_##_layer##_pdelab_cg_##_p##_double_##_r##_##_rC =                                  \
//      Dune::GDT::bindings::SystemAssembler<Dune::GDT::CgSpaceProvider<_GRID,                                         \
//                                                                      Dune::XT::Grid::Layers::_layer,                \
//                                                                      Dune::GDT::Backends::pdelab,         \
//                                                                      _p,                                            \
//                                                                      double,                                        \
//                                                                      _r,                                            \
//                                                                      _rC>>::bind(_m);                               \
//  Dune::GDT::bindings::                                                                                              \
//      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                   \
//                               typename Dune::XT::Grid::Layer<_GRID,                                                 \
//                                                              Dune::XT::Grid::Layers::_layer,                        \
//                                                              Dune::XT::Grid::Backends::view,                        \
//                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,      \
//                           _GRID>::addbind(system_assembler_##_GRID##_##_layer##_pdelab_cg_##_p##_double_##_r##_##_rC)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALBERTA(_m, _layer, _p)                                            \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG(_m, ALBERTA_2D, _layer, _p, 1, 1)
//#else
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALBERTA(_m, _p, _layer)
//#endif

//#if HAVE_DUNE_ALUGRID
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALU(_m, _layer, _p)                                                \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1)
//#else
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALU(_m, _layer, _p)
//#endif

////#if HAVE_DUNE_UGGRID || HAVE_UG // <- does not work
////#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_UG(_m, _layer, _p)
////  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG(_m, UG_2D, _layer, _p, 1, 1)
////#else
////#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_UG(_m, _layer, _p)
////#endif

//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_YASP(_m, _layer, _p)                                               \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1)

//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALL(_m, _layer, _p)                                                \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALBERTA(_m, _layer, _p);                                                 \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALU(_m, _layer, _p);                                                     \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_YASP(_m, _layer, _p)
////_DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_UG(_m, _layer, _p); // <- does not work

//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(_m)                                                                   \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALL(_m, leaf, 1);                                                        \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB_CG_ALL(_m, level, 1);

//#else // HAVE_DUNE_PDELAB
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(_m)
//#endif

// * gdt backend
//   - fv spaces

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV(_m, _GRID, _layer, _r, _rC)                                             \
  auto system_assembler_##_GRID##_##_layer##_gdt_fv_double_##_r##_##_rC =                                              \
      Dune::GDT::bindings::SystemAssembler<Dune::GDT::FvSpaceProvider<_GRID,                                           \
                                                                      Dune::XT::Grid::Layers::_layer,                  \
                                                                      Dune::GDT::Backends::gdt,                        \
                                                                      double,                                          \
                                                                      _r,                                              \
                                                                      _rC>>::bind(_m);                                 \
  Dune::GDT::bindings::                                                                                                \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::view,                          \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::addbind(system_assembler_##_GRID##_##_layer##_gdt_fv_double_##_r##_##_rC)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALBERTA(_m, _layer)                                                   \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV(_m, ALBERTA_2D, _layer, 1, 1)
//#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALBERTA(_m, _layer)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALU(_m, _layer)                                                         \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, 1, 1)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALU(_m, _layer)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_UG(_m, _layer)                                                        \
//  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV(_m, UG_2D, _layer, 1, 1)
//#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_UG(_m, _layer)
//#endif

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_YASP(_m, _layer)                                                        \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, 1, 1)

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALL(_m, _layer)                                                         \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALBERTA(_m, _layer);                                                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALU(_m, _layer);                                                              \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_UG(_m, _layer);                                                               \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_YASP(_m, _layer)

#define _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(_m)                                                                        \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALL(_m, leaf);                                                                \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT_FV_ALL(_m, level)

#define DUNE_GDT_ASSEMBLER_SYSTEM_BIND(_m)                                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(_m);                                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(_m);                                                                             \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
