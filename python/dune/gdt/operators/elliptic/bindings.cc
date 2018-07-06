// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>

#include <python/dune/gdt/operators/elliptic/bindings.hh>


#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND(_m, _G, _gl, _gl_backend, _dt, _la, _s_backend, _s_type, _s_layer, _s_p)      \
  Dune::GDT::bindings::EllipticMatrixOperator<_G,                                                                      \
                                              Dune::XT::Grid::Layers::_gl,                                             \
                                              Dune::XT::Grid::Backends::_gl_backend,                                   \
                                              _dt,                                                                     \
                                              Dune::XT::LA::Backends::_la,                                             \
                                              Dune::GDT::Backends::_s_backend,                                         \
                                              Dune::GDT::SpaceType::_s_type,                                           \
                                              Dune::XT::Grid::Layers::_s_layer,                                        \
                                              _s_p>::bind(_m)

PYBIND11_MODULE(__operators_elliptic, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune;
  using G = GDT_BINDINGS_GRID;
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.operators.elliptic");
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, view, true, istl_sparse, gdt, dg, leaf, 1);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, view, true, istl_sparse, gdt, dg, leaf, 2);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, view, true, istl_sparse, gdt, dg, leaf, 3);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, dd_subdomain, view, true, istl_sparse, gdt, dg, dd_subdomain, 1);
}
