// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef PYTHON_DUNE_GDT_PROJECTIONS_BINDINGS_HH
#define PYTHON_DUNE_GDT_PROJECTIONS_BINDINGS_HH
//#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces/bindings.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/projections.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class SP, class V>
class project
{
  typedef typename SP::type SpaceType;
  static_assert(is_space<SpaceType>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");
  typedef typename XT::Functions::LocalizableFunctionInterface<typename SpaceType::EntityType,
                                                               typename SpaceType::DomainFieldType,
                                                               SpaceType::dimDomain,
                                                               typename SpaceType::RangeFieldType,
                                                               SpaceType::dimRange,
                                                               SpaceType::dimRangeCols>
      SourceType;
  typedef DiscreteFunction<SpaceType, V> RangeType;

public:
  static void bind(pybind11::module& m)
  {
    using namespace pybind11::literals;

    m.def("project",
          [](const SourceType& source, RangeType& range, const size_t over_integrate) {
            GDT::project(source, range, over_integrate);
          },
          "source"_a,
          "range"_a,
          "over_integrate"_a = 0);
  } // ... bind(...)
}; // class project


// begin: this is what we need for the .so

#define _DUNE_GDT_PROJECTIONS_BIND(_m, _GRID, _g_layer, _s_backend, _s_type, _p, _la)                                  \
  Dune::GDT::bindings::project<Dune::GDT::SpaceProvider<_GRID,                                                         \
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
#define _DUNE_GDT_PROJECTIONS_BIND_ALBERTA(_m, _g_layer, _s_backend, _s_type, _p, _la)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_PROJECTIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la)                                     \
  _DUNE_GDT_PROJECTIONS_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, _g_layer, _s_backend, _s_type, _p, _la)
#else
#define _DUNE_GDT_PROJECTIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//  ...
//#else
#define _DUNE_GDT_PROJECTIONS_BIND_UG(_m, _g_layer, _s_backend, _s_type, _p, _la)
//#endif

#define _DUNE_GDT_PROJECTIONS_BIND_YASP(_m, _g_layer, _s_backend, _s_type, _p, _la)                                    \
  _DUNE_GDT_PROJECTIONS_BIND(_m, YASP_1D_EQUIDISTANT_OFFSET, _g_layer, _s_backend, _s_type, _p, _la);                  \
  _DUNE_GDT_PROJECTIONS_BIND(_m, YASP_2D_EQUIDISTANT_OFFSET, _g_layer, _s_backend, _s_type, _p, _la)

#define _DUNE_GDT_PROJECTIONS_BIND_ALL_GRIDS(_m, _g_layer, _s_backend, _s_type, _p, _la)                               \
  _DUNE_GDT_PROJECTIONS_BIND_ALBERTA(_m, _g_layer, _s_backend, _s_type, _p, _la);                                      \
  _DUNE_GDT_PROJECTIONS_BIND_ALU(_m, _g_layer, _s_backend, _s_type, _p, _la);                                          \
  _DUNE_GDT_PROJECTIONS_BIND_UG(_m, _g_layer, _s_backend, _s_type, _p, _la);                                           \
  _DUNE_GDT_PROJECTIONS_BIND_YASP(_m, _g_layer, _s_backend, _s_type, _p, _la)

// for each space backend

#define _DUNE_GDT_PROJECTIONS_BIND_DEFAULT(_m, _la)                                                                    \
  _DUNE_GDT_PROJECTIONS_BIND_ALL_GRIDS(_m, leaf, gdt, fv, 0, _la);                                                     \
  _DUNE_GDT_PROJECTIONS_BIND_ALL_GRIDS(_m, dd_subdomain, gdt, cg, 1, _la);                                             \
  _DUNE_GDT_PROJECTIONS_BIND_ALL_GRIDS(_m, level, gdt, fv, 0, _la);                                                    \
  _DUNE_GDT_PROJECTIONS_BIND_ALL_GRIDS(_m, dd_subdomain, gdt, dg, 1, _la)


#define _DUNE_GDT_PROJECTIONS_BIND_ALL_SPACES(_m, _la) _DUNE_GDT_PROJECTIONS_BIND_DEFAULT(_m, _la);

// for each la backend

//#define _DUNE_GDT_PROJECTIONS_BIND_COMMON(_m) _DUNE_GDT_PROJECTIONS_BIND_ALL_SPACES(_m, common_dense)
#define _DUNE_GDT_PROJECTIONS_BIND_COMMON(_m)

//#if HAVE_EIGEN
//#define _DUNE_GDT_PROJECTIONS_BIND_EIGEN(_m) _DUNE_GDT_PROJECTIONS_BIND_ALL_SPACES(_m, eigen_dense)
//#else
#define _DUNE_GDT_PROJECTIONS_BIND_EIGEN(_m)
//#endif

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_PROJECTIONS_BIND_ISTL(_m) _DUNE_GDT_PROJECTIONS_BIND_ALL_SPACES(_m, istl_dense)
#else
#define _DUNE_GDT_PROJECTIONS_BIND_ISTL(_m)
#endif

#define DUNE_GDT_PROJECTIONS_BIND(_m)                                                                                  \
  _DUNE_GDT_PROJECTIONS_BIND_COMMON(_m);                                                                               \
  _DUNE_GDT_PROJECTIONS_BIND_EIGEN(_m);                                                                                \
  _DUNE_GDT_PROJECTIONS_BIND_ISTL(_m)

// end: this is what we need for the .so


} // namespace bindings
} // namespace GDT
} // namespace Dune

//#endif // HAVE_DUNE_PYBINDXI
#endif // PYTHON_DUNE_GDT_PROJECTIONS_BINDINGS_HH
