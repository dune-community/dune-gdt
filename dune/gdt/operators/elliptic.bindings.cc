// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>

#include <dune/gdt/operators/elliptic.bindings.hh>


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


PYBIND11_PLUGIN(__operators_elliptic)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune;

  py::module m("__operators_elliptic", "dune-gdt: EllipticMatrixOperator");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.operators.elliptic");

  using G = ALU_2D_SIMPLEX_CONFORMING;

  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, part, true, istl_sparse, fem, dg, leaf, 1);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, part, true, istl_sparse, fem, dg, leaf, 2);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, leaf, part, true, istl_sparse, fem, dg, leaf, 3);
  DUNE_GDT_OPERATORS_ELLIPTIC_BIND(m, G, dd_subdomain, part, true, istl_sparse, fem, dg, dd_subdomain, 1);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
