// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>

#include <dune/gdt/operators/elliptic-ipdg.bindings.hh>


#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(                                                                         \
    _m, _G, _gl, _gl_backend, _dt, _ipdg, _la, _s_backend, _s_type, _s_layer, _s_p)                                    \
  Dune::GDT::bindings::EllipticIpdgMatrixOperator<_G,                                                                  \
                                                  Dune::XT::Grid::Layers::_gl,                                         \
                                                  Dune::XT::Grid::Backends::_gl_backend,                               \
                                                  _dt,                                                                 \
                                                  Dune::GDT::LocalEllipticIpdgIntegrands::Method::_ipdg,               \
                                                  Dune::XT::LA::Backends::_la,                                         \
                                                  Dune::GDT::Backends::_s_backend,                                     \
                                                  Dune::GDT::SpaceType::_s_type,                                       \
                                                  Dune::XT::Grid::Layers::_s_layer,                                    \
                                                  _s_p>::bind(_m)

PYBIND11_MODULE(__local_elliptic_ipdg_operators, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using G = ALU_2D_SIMPLEX_CONFORMING;
  add_initialization(m, "dune.gdt.operators.elliptic.ipdg");
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(m, G, leaf, part, true, swipdg_affine_factor, istl_sparse, gdt, dg, leaf, 1);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(m, G, leaf, part, true, swipdg_affine_factor, istl_sparse, gdt, dg, leaf, 2);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(m, G, leaf, part, true, swipdg_affine_factor, istl_sparse, gdt, dg, leaf, 3);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(
      m, G, dd_subdomain, part, true, swipdg_affine_factor, istl_sparse, gdt, dg, dd_subdomain, 1);
}

#endif // HAVE_DUNE_PYBINDXI
