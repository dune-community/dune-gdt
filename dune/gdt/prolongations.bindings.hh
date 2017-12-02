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

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.hh>
#include <dune/gdt/playground/spaces/block.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "prolongations.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class SP, class V>
class prolong
{
  typedef typename SP::type SpaceType;
  static_assert(is_space<SpaceType>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");
  typedef ConstDiscreteFunction<SpaceType, V> SourceType;
  typedef DiscreteFunction<SpaceType, V> RangeType;

public:
  static void bind(pybind11::module& m)
  {
    using namespace pybind11::literals;

    m.def("prolong",
          [](const SourceType& source, RangeType& range, const size_t over_integrate) {
            GDT::prolong(source, range, over_integrate);
          },
          "source"_a,
          "range"_a,
          "over_integrate"_a = 0);
  } // ... bind(...)
}; // class prolong


// begin: this is what we need for the .so

#define _DUNE_GDT_PROLONGATIONS_BIND(_m, _GRID, _g_layer, _s_backend, _s_type, _p, _la)                                \
  Dune::GDT::bindings::prolong<Dune::GDT::SpaceProvider<_GRID,                                                         \
                                                        Dune::XT::Grid::Layers::_g_layer,                              \
                                                        Dune::GDT::SpaceType::_s_type,                                 \
                                                        Dune::GDT::Backends::_s_backend,                               \
                                                        _p,                                                            \
                                                        double,                                                        \
                                                        1,                                                             \
                                                        1>,                                                            \
                               typename Dune::XT::LA::Container<double,                                                \
                                                                Dune::XT::LA::Backends::_la>::VectorType>::bind(_m)

// for each grid

//#if HAVE_ALBERTA
//  ...
//#else
#define _DUNE_GDT_PROLONGATIONS_BIND_ALBERTA(_m, _g_layer, _s_backend, _s_type, _p, _la)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_PROLONGATIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la)                                   \
  _DUNE_GDT_PROLONGATIONS_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, _g_layer, _s_backend, _s_type, _p, _la)
#else
#define _DUNE_GDT_PROLONGATIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//  ...
//#else
#define _DUNE_GDT_PROLONGATIONS_BIND_UG(_m, _g_layer, _s_backend, _s_type, _p, _la)
//#endif

#define _DUNE_GDT_PROLONGATIONS_BIND_YASP(_m, _g_layer, _s_backend, _s_type, _p, _la)
//  _DUNE_GDT_PROLONGATIONS_BIND(_m, YASP_1D_EQUIDISTANT_OFFSET, _g_layer, _s_backend, _s_type, _p, _la);              \
//  _DUNE_GDT_PROLONGATIONS_BIND(_m, YASP_2D_EQUIDISTANT_OFFSET, _g_layer, _s_backend, _s_type, _p, _la)

#define _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, _g_layer, _s_backend, _s_type, _p, _la)                             \
  _DUNE_GDT_PROLONGATIONS_BIND_ALBERTA(_m, _g_layer, _s_backend, _s_type, _p, _la);                                    \
  _DUNE_GDT_PROLONGATIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la);                                        \
  _DUNE_GDT_PROLONGATIONS_BIND_UG(_m, _g_layer, _s_backend, _s_type, _p, _la);                                         \
  _DUNE_GDT_PROLONGATIONS_BIND_YASP(_m, _g_layer, _s_backend, _s_type, _p, _la)

// for each space backend

#define _DUNE_GDT_PROLONGATIONS_BIND_DEFAULT(_m, _la)                                                                  \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, leaf, gdt, fv, 0, _la);                                                   \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, level, gdt, fv, 0, _la)

#if HAVE_DUNE_FEM
#define _DUNE_GDT_PROLONGATIONS_BIND_FEM(_m, _la)                                                                      \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, leaf, fem, cg, 1, _la);                                                   \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, level, fem, cg, 1, _la);                                                  \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, dd_subdomain, fem, block_cg, 1, _la);                                     \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, leaf, fem, dg, 1, _la);                                                   \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, level, fem, dg, 1, _la);                                                  \
  _DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, dd_subdomain, fem, block_dg, 1, _la)
// These do not work with alugrid:
//_DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, dd_subdomain, fem, cg, 1, _la);                                           \
//_DUNE_GDT_PROLONGATIONS_BIND_ALL_GRIDS(_m, dd_subdomain, fem, dg, 1, _la);                                           \

#else
#define _DUNE_GDT_PROLONGATIONS_BIND_FEM(_m, _la)
#endif

//#if HAVE_DUNE_FUNCTIONS
//  ...
//#else
#define _DUNE_GDT_PROLONGATIONS_BIND_FUNCTIONS(_m, _la)
//#endif

//#if HAVE_DUNE_PDELAB
//  ...
//#else
#define _DUNE_GDT_PROLONGATIONS_BIND_PDELAB(_m, _la)
//#endif

#define _DUNE_GDT_PROLONGATIONS_BIND_ALL_SPACES(_m, _la)                                                               \
  _DUNE_GDT_PROLONGATIONS_BIND_DEFAULT(_m, _la);                                                                       \
  _DUNE_GDT_PROLONGATIONS_BIND_FEM(_m, _la);                                                                           \
  _DUNE_GDT_PROLONGATIONS_BIND_FUNCTIONS(_m, _la);                                                                     \
  _DUNE_GDT_PROLONGATIONS_BIND_PDELAB(_m, _la)

// for each la backend

//#define _DUNE_GDT_PROLONGATIONS_BIND_COMMON(_m) _DUNE_GDT_PROLONGATIONS_BIND_ALL_SPACES(_m, common_dense)
#define _DUNE_GDT_PROLONGATIONS_BIND_COMMON(_m)

//#if HAVE_EIGEN
//#define _DUNE_GDT_PROLONGATIONS_BIND_EIGEN(_m) _DUNE_GDT_PROLONGATIONS_BIND_ALL_SPACES(_m, eigen_dense)
//#else
#define _DUNE_GDT_PROLONGATIONS_BIND_EIGEN(_m)
//#endif

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_PROLONGATIONS_BIND_ISTL(_m) _DUNE_GDT_PROLONGATIONS_BIND_ALL_SPACES(_m, istl_dense)
#else
#define _DUNE_GDT_PROLONGATIONS_BIND_ISTL(_m)
#endif

#define DUNE_GDT_PROLONGATIONS_BIND(_m) _DUNE_GDT_PROLONGATIONS_BIND_ISTL(_m)
//  _DUNE_GDT_PROLONGATIONS_BIND_COMMON(_m);                                                                           \
//  _DUNE_GDT_PROLONGATIONS_BIND_EIGEN(_m);                                                                            \

// end: this is what we need for the .so


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_PROLONGATIONS_BINDINGS_HH
