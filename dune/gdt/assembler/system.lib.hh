// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_LIB_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_LIB_HH

#if DUNE_XT_WITH_PYTHON_BINDINGS

#include <dune/xt/la/container.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/spaces.hh>

#include "system.hh"


// everything related to dd subdomain
#if HAVE_DUNE_FEM
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN(                                                            \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, _la)                              \
  _pre Dune::GDT::SystemAssembler<typename Dune::GDT::SpaceProvider<_G,                                                \
                                                                    Dune::XT::Grid::Layers::_s_grid_layer,             \
                                                                    Dune::GDT::SpaceType::_s_type,                     \
                                                                    Dune::GDT::Backends::_s_backend,                   \
                                                                    _p,                                                \
                                                                    double,                                            \
                                                                    _r,                                                \
                                                                    _rC>::type,                                        \
                                  typename Dune::XT::Grid::Layer<_G,                                                   \
                                                                 Dune::XT::Grid::Layers::_g_layer,                     \
                                                                 Dune::XT::Grid::Backends::_g_backend,                 \
                                                                 Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>&        \
  Dune::GDT::SystemAssembler<typename Dune::GDT::SpaceProvider<_G,                                                     \
                                                               Dune::XT::Grid::Layers::_s_grid_layer,                  \
                                                               Dune::GDT::SpaceType::_s_type,                          \
                                                               Dune::GDT::Backends::_s_backend,                        \
                                                               _p,                                                     \
                                                               double,                                                 \
                                                               _r,                                                     \
                                                               _rC>::type,                                             \
                             typename Dune::XT::Grid::Layer<_G,                                                        \
                                                            Dune::XT::Grid::Layers::_g_layer,                          \
                                                            Dune::XT::Grid::Backends::_g_backend,                      \
                                                            Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>::            \
      append<typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType::Traits, double>(       \
          const Dune::GDT::LocalVolumeTwoFormAssembler<                                                                \
              typename Dune::GDT::SpaceProvider<_G,                                                                    \
                                                Dune::XT::Grid::Layers::_s_grid_layer,                                 \
                                                Dune::GDT::SpaceType::_s_type,                                         \
                                                Dune::GDT::Backends::_s_backend,                                       \
                                                _p,                                                                    \
                                                double,                                                                \
                                                _r,                                                                    \
                                                _rC>::type,                                                            \
              typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType,                       \
              typename Dune::GDT::SpaceProvider<_G,                                                                    \
                                                Dune::XT::Grid::Layers::_s_grid_layer,                                 \
                                                Dune::GDT::SpaceType::_s_type,                                         \
                                                Dune::GDT::Backends::_s_backend,                                       \
                                                _p,                                                                    \
                                                double,                                                                \
                                                _r,                                                                    \
                                                _rC>::type>&,                                                          \
          Dune::XT::LA::MatrixInterface<                                                                               \
              typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType::Traits,               \
              double>&,                                                                                                \
          const XT::Grid::ApplyOn::WhichEntity<                                                                        \
              typename Dune::XT::Grid::Layer<_G,                                                                       \
                                             Dune::XT::Grid::Layers::_g_layer,                                         \
                                             Dune::XT::Grid::Backends::_g_backend,                                     \
                                             Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>*)

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN_ISTL(                                                       \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)                                   \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN(                                                                  \
      _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, istl_sparse)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN_ISTL(                                                       \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)
#endif

#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                   \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)                                   \
  _pre class Dune::GDT::SystemAssembler<typename Dune::GDT::SpaceProvider<_G,                                          \
                                                                          Dune::XT::Grid::Layers::_s_grid_layer,       \
                                                                          Dune::GDT::SpaceType::_s_type,               \
                                                                          Dune::GDT::Backends::_s_backend,             \
                                                                          _p,                                          \
                                                                          double,                                      \
                                                                          _r,                                          \
                                                                          _rC>::type,                                  \
                                        typename Dune::XT::Grid::Layer<_G,                                             \
                                                                       Dune::XT::Grid::Layers::_g_layer,               \
                                                                       Dune::XT::Grid::Backends::_g_backend,           \
                                                                       Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>;  \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN_ISTL(                                                             \
      _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                   \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)
#endif

// everything not related to dd subdomain
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND(                                                                         \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, _la)                              \
  _pre Dune::GDT::SystemAssembler<typename Dune::GDT::SpaceProvider<_G,                                                \
                                                                    Dune::XT::Grid::Layers::_s_grid_layer,             \
                                                                    Dune::GDT::SpaceType::_s_type,                     \
                                                                    Dune::GDT::Backends::_s_backend,                   \
                                                                    _p,                                                \
                                                                    double,                                            \
                                                                    _r,                                                \
                                                                    _rC>::type,                                        \
                                  typename Dune::XT::Grid::Layer<_G,                                                   \
                                                                 Dune::XT::Grid::Layers::_g_layer,                     \
                                                                 Dune::XT::Grid::Backends::_g_backend>::type>&         \
  Dune::GDT::SystemAssembler<                                                                                          \
      typename Dune::GDT::SpaceProvider<_G,                                                                            \
                                        Dune::XT::Grid::Layers::_s_grid_layer,                                         \
                                        Dune::GDT::SpaceType::_s_type,                                                 \
                                        Dune::GDT::Backends::_s_backend,                                               \
                                        _p,                                                                            \
                                        double,                                                                        \
                                        _r,                                                                            \
                                        _rC>::type,                                                                    \
      typename Dune::XT::Grid::Layer<_G, Dune::XT::Grid::Layers::_g_layer, Dune::XT::Grid::Backends::_g_backend>::     \
          type>::append<typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType::Traits,     \
                        double>(                                                                                       \
      const Dune::GDT::LocalVolumeTwoFormAssembler<                                                                    \
          typename Dune::GDT::SpaceProvider<_G,                                                                        \
                                            Dune::XT::Grid::Layers::_s_grid_layer,                                     \
                                            Dune::GDT::SpaceType::_s_type,                                             \
                                            Dune::GDT::Backends::_s_backend,                                           \
                                            _p,                                                                        \
                                            double,                                                                    \
                                            _r,                                                                        \
                                            _rC>::type,                                                                \
          typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType,                           \
          typename Dune::GDT::SpaceProvider<_G,                                                                        \
                                            Dune::XT::Grid::Layers::_s_grid_layer,                                     \
                                            Dune::GDT::SpaceType::_s_type,                                             \
                                            Dune::GDT::Backends::_s_backend,                                           \
                                            _p,                                                                        \
                                            double,                                                                    \
                                            _r,                                                                        \
                                            _rC>::type>&,                                                              \
      Dune::XT::LA::MatrixInterface<                                                                                   \
          typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::MatrixType::Traits,                   \
          double>&,                                                                                                    \
      const XT::Grid::ApplyOn::WhichEntity<                                                                            \
          typename Dune::XT::Grid::Layer<_G, Dune::XT::Grid::Layers::_g_layer, Dune::XT::Grid::Backends::_g_backend>:: \
              type>*)

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_ISTL(                                                                    \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)                                   \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND(                                                                               \
      _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, istl_sparse)
#else
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_ISTL(                                                                    \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)
#endif

#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB(                                                                                \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC)                                   \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_ISTL(                                                                          \
      _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC);                                \
  _pre class Dune::GDT::SystemAssembler<typename Dune::GDT::SpaceProvider<_G,                                          \
                                                                          Dune::XT::Grid::Layers::_s_grid_layer,       \
                                                                          Dune::GDT::SpaceType::_s_type,               \
                                                                          Dune::GDT::Backends::_s_backend,             \
                                                                          _p,                                          \
                                                                          double,                                      \
                                                                          _r,                                          \
                                                                          _rC>::type,                                  \
                                        typename Dune::XT::Grid::Layer<_G,                                             \
                                                                       Dune::XT::Grid::Layers::_g_layer,               \
                                                                       Dune::XT::Grid::Backends::_g_backend>::type>

#define DUNE_GDT_ASSEMBLER_SYSTEM_LIB(_pre, _G, _s_type, _s_backend, _p, _r, _rC)                                      \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB(_pre, _G, leaf, view, _s_type, _s_backend, leaf, _p, _r, _rC);                        \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB(_pre, _G, level, view, _s_type, _s_backend, level, _p, _r, _rC)

#define DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(_pre, _G, _s_type, _s_backend, _p, _r, _rC)                         \
  DUNE_GDT_ASSEMBLER_SYSTEM_LIB(_pre, _G, _s_type, _s_backend, _p, _r, _rC);                                           \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, adaptive_leaf, part, _s_type, _s_backend, adaptive_leaf, _p, _r, _rC);                                 \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(_pre, _G, leaf, part, _s_type, _s_backend, leaf, _p, _r, _rC);           \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(_pre, _G, level, part, _s_type, _s_backend, level, _p, _r, _rC);         \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain, part, _s_type, _s_backend, dd_subdomain, _p, _r, _rC);                                   \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain_boundary, part, _s_type, _s_backend, dd_subdomain, _p, _r, _rC);                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain_coupling, part, _s_type, _s_backend, dd_subdomain, _p, _r, _rC)


#if HAVE_DUNE_FEM
#if HAVE_DUNE_ALUGRID
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, ALU_2D_SIMPLEX_CONFORMING, cg, fem, 1, 1, 1);
#endif
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_1D_EQUIDISTANT_OFFSET, cg, fem, 1, 1, 1);
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_2D_EQUIDISTANT_OFFSET, cg, fem, 1, 1, 1);
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_3D_EQUIDISTANT_OFFSET, cg, fem, 1, 1, 1);
#endif


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
#endif // DUNE_GDT_ASSEMBLER_SYSTEM_LIB_HH
