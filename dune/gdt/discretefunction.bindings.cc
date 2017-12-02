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

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>

#include <dune/gdt/discretefunction/default.bindings.hh>


#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m, _G, _g_layer, _s_backend, _s_type, _p, _r, _la)                     \
  Dune::GDT::bindings::                                                                                                \
      DiscreteFunction<Dune::GDT::SpaceProvider<_G,                                                                    \
                                                Dune::XT::Grid::Layers::_g_layer,                                      \
                                                Dune::GDT::SpaceType::_s_type,                                         \
                                                Dune::GDT::Backends::_s_backend,                                       \
                                                _p,                                                                    \
                                                double,                                                                \
                                                _r,                                                                    \
                                                1>,                                                                    \
                       typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::VectorType>::bind(_m)


PYBIND11_PLUGIN(__discretefunction)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__discretefunction", "dune-gdt: discrete functions");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.discretefunction");

  using G = ALU_2D_SIMPLEX_CONFORMING;

  DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(m, G, leaf, fem, dg, 1, 1, istl_dense);
  DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(m, G, leaf, fem, dg, 2, 1, istl_dense);
  DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(m, G, leaf, fem, dg, 3, 1, istl_dense);
  DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(m, G, dd_subdomain, fem, block_dg, 1, 1, istl_dense);
  DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(m, G, leaf, pdelab, rt, 0, 2, istl_dense);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
