// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
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
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/spaces.hh>

#include "system.hh"


// everything related to dd subdomain
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND_DD_SUBDOMAIN(                                                                                                           \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, _la)                                                                             \
  typedef typename Dune::GDT::SpaceProvider<_G,                                                                                                                       \
                                            Dune::XT::Grid::Layers::_s_grid_layer,                                                                                    \
                                            Dune::GDT::SpaceType::_s_type,                                                                                            \
                                            Dune::GDT::Backends::_s_backend,                                                                                          \
                                            _p,                                                                                                                       \
                                            double,                                                                                                                   \
                                            _r,                                                                                                                       \
                                            _rC>::type                                                                                                                \
      _DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                        \
  typedef typename _DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                          \
      BaseFunctionSetType                                                                                                                                             \
          _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                          \
  typedef typename Dune::XT::Grid::Layer<_G,                                                                                                                          \
                                         Dune::XT::Grid::Layers::_g_layer,                                                                                            \
                                         Dune::XT::Grid::Backends::_g_backend,                                                                                        \
                                         Dune::XT::Grid::DD::SubdomainGrid<_G>>::type                                                                                 \
      _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                    \
  typedef Dune::XT::Functions::                                                                                                                                       \
      LocalizableFunctionInterface<Dune::XT::Grid::                                                                                                                   \
                                       extract_entity_t<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>, \
                                   typename _G::ctype,                                                                                                                \
                                   _G::dimension,                                                                                                                     \
                                   double,                                                                                                                            \
                                   1,                                                                                                                                 \
                                   1>                                                                                                                                 \
          _DUNE_GDT_ASLADS_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                           \
  typedef Dune::XT::Functions::                                                                                                                                       \
      LocalfunctionInterface<Dune::XT::Grid::                                                                                                                         \
                                 extract_entity_t<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>,       \
                             typename _G::ctype,                                                                                                                      \
                             _G::dimension,                                                                                                                           \
                             double,                                                                                                                                  \
                             1,                                                                                                                                       \
                             1>                                                                                                                                       \
          _DUNE_GDT_ASLADS_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                      \
  typedef Dune::XT::Grid::                                                                                                                                            \
      extract_intersection_t<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>                             \
          _DUNE_GDT_ASLADS_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                             \
  typedef typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::istl_sparse>::MatrixType                                                                   \
      _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                         \
          const Dune::GDT::                                                                                                                                           \
              LocalVolumeTwoFormInterface<_DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,          \
                                          _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,          \
                                          double>&,                                                                                                                   \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          const XT::Grid::ApplyOn::                                                                                                                                   \
              WhichEntity<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                             \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                         \
          const Dune::GDT::                                                                                                                                           \
              LocalCouplingTwoFormInterface<_DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                                 \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          const Dune::XT::Grid::ApplyOn::                                                                                                                             \
              WhichIntersection<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                         \
          const Dune::GDT::                                                                                                                                           \
              LocalCouplingTwoFormInterface<_DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                                 \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          const Dune::XT::Grid::ApplyOn::                                                                                                                             \
              WhichIntersection<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                         \
          const Dune::GDT::                                                                                                                                           \
              LocalBoundaryTwoFormInterface<_DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLADS_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLADS_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                                 \
          Dune::XT::LA::MatrixInterface<                                                                                                                              \
              typename _DUNE_GDT_ASLADS_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                             \
              double>&,                                                                                                                                               \
          const Dune::XT::Grid::ApplyOn::                                                                                                                             \
              WhichIntersection<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLADS_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                         \
          const Dune::GDT::                                                                                                                                           \
              LocalVolumeTwoFormInterface<_DUNE_GDT_ASLADS_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,      \
                                          _DUNE_GDT_ASLADS_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,      \
                                          double>&,                                                                                                                   \
          const _DUNE_GDT_ASLADS_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la&,                                    \
          const _DUNE_GDT_ASLADS_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la&,                                    \
          double&,                                                                                                                                                    \
          const XT::Grid::ApplyOn::                                                                                                                                   \
              WhichEntity<_DUNE_GDT_ASLADS_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*)

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

// everything not related to dd subdomain
#define _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_APPEND(                                                                                                                      \
    _pre, _G, _g_layer, _g_backend, _s_type, _s_backend, _s_grid_layer, _p, _r, _rC, _la)                                                                           \
  typedef typename Dune::GDT::SpaceProvider<_G,                                                                                                                     \
                                            Dune::XT::Grid::Layers::_s_grid_layer,                                                                                  \
                                            Dune::GDT::SpaceType::_s_type,                                                                                          \
                                            Dune::GDT::Backends::_s_backend,                                                                                        \
                                            _p,                                                                                                                     \
                                            double,                                                                                                                 \
                                            _r,                                                                                                                     \
                                            _rC>::type                                                                                                              \
      _DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                        \
  typedef typename _DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                          \
      BaseFunctionSetType                                                                                                                                           \
          _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                          \
  typedef                                                                                                                                                           \
      typename Dune::XT::Grid::Layer<_G, Dune::XT::Grid::Layers::_g_layer, Dune::XT::Grid::Backends::_g_backend>::type                                              \
          _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                \
  typedef Dune::XT::Functions::                                                                                                                                     \
      LocalizableFunctionInterface<Dune::XT::Grid::                                                                                                                 \
                                       extract_entity_t<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>, \
                                   typename _G::ctype,                                                                                                              \
                                   _G::dimension,                                                                                                                   \
                                   double,                                                                                                                          \
                                   1,                                                                                                                               \
                                   1>                                                                                                                               \
          _DUNE_GDT_ASLA_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                           \
  typedef Dune::XT::Functions::                                                                                                                                     \
      LocalfunctionInterface<Dune::XT::Grid::                                                                                                                       \
                                 extract_entity_t<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>,       \
                             typename _G::ctype,                                                                                                                    \
                             _G::dimension,                                                                                                                         \
                             double,                                                                                                                                \
                             1,                                                                                                                                     \
                             1>                                                                                                                                     \
          _DUNE_GDT_ASLA_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                      \
  typedef Dune::XT::Grid::                                                                                                                                          \
      extract_intersection_t<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>                             \
          _DUNE_GDT_ASLA_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                             \
  typedef typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::istl_sparse>::MatrixType                                                                 \
      _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la;                                                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                       \
          const Dune::GDT::                                                                                                                                         \
              LocalVolumeTwoFormInterface<_DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,          \
                                          _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,          \
                                          double>&,                                                                                                                 \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          const XT::Grid::ApplyOn::                                                                                                                                 \
              WhichEntity<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                             \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                       \
          const Dune::GDT::                                                                                                                                         \
              LocalCouplingTwoFormInterface<_DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          const Dune::XT::Grid::ApplyOn::                                                                                                                           \
              WhichIntersection<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                       \
          const Dune::GDT::                                                                                                                                         \
              LocalCouplingTwoFormInterface<_DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          const Dune::XT::Grid::ApplyOn::                                                                                                                           \
              WhichIntersection<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                       \
          const Dune::GDT::                                                                                                                                         \
              LocalBoundaryTwoFormInterface<_DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            _DUNE_GDT_ASLA_Intersection_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,           \
                                            _DUNE_GDT_ASLA_BaseFunctionSet_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,        \
                                            double>&,                                                                                                               \
          Dune::XT::LA::MatrixInterface<                                                                                                                            \
              typename _DUNE_GDT_ASLA_Matrix_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la::                                     \
                  Traits,                                                                                                                                           \
              double>&,                                                                                                                                             \
          const Dune::XT::Grid::ApplyOn::                                                                                                                           \
              WhichIntersection<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*);                       \
  _pre Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                            \
                                  _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>&                       \
  Dune::GDT::SystemAssembler<_DUNE_GDT_ASLA_Space_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,                                 \
                             _DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>::                           \
      append(                                                                                                                                                       \
          const Dune::GDT::                                                                                                                                         \
              LocalVolumeTwoFormInterface<_DUNE_GDT_ASLA_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,      \
                                          _DUNE_GDT_ASLA_LocalScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la,      \
                                          double>&,                                                                                                                 \
          const _DUNE_GDT_ASLA_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la&,                                    \
          const _DUNE_GDT_ASLA_ScalarFunction_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la&,                                    \
          double&,                                                                                                                                                  \
          const XT::Grid::ApplyOn::                                                                                                                                 \
              WhichEntity<_DUNE_GDT_ASLA_GridLayer_##_G##_g_layer##_g_backend##_s_type##_s_backend##_s_grid_layer##_p##_r##_rC##_la>*)

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
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(_pre, _G, leaf, view, _s_type, _s_backend, leaf, _p, _r, _rC);           \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(_pre, _G, level, view, _s_type, _s_backend, level, _p, _r, _rC);         \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain, view, _s_type, _s_backend, dd_subdomain, _p, _r, _rC);                                   \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain_boundary, view, _s_type, _s_backend, dd_subdomain, _p, _r, _rC);                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(                                                                         \
      _pre, _G, dd_subdomain_coupling, view, _s_type, _s_backend, dd_subdomain, _p, _r, _rC)
//DUNE_GDT_ASSEMBLER_SYSTEM_LIB(_pre, _G, _s_type, _s_backend, _p, _r, _rC);                                           \


#if HAVE_DUNE_ALUGRID
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, ALU_2D_SIMPLEX_CONFORMING, cg, gdt, 1, 1, 1);
#endif
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_1D_EQUIDISTANT_OFFSET, cg, gdt, 1, 1, 1);
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_2D_EQUIDISTANT_OFFSET, cg, gdt, 1, 1, 1);
DUNE_GDT_ASSEMBLER_SYSTEM_LIB_DD_SUBDOMAIN(extern template, YASP_3D_EQUIDISTANT_OFFSET, cg, gdt, 1, 1, 1);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
#endif // DUNE_GDT_ASSEMBLER_SYSTEM_LIB_HH
