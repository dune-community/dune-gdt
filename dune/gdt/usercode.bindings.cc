// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

/**
  * This file is intended as a quick starting point to add Python bindings for user code.
  **/


#include "config.h"

#include <csignal>

#if HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/bindings.hh>

#include <dune/xt/la/container.bindings.hh>
#include <dune/xt/la/container/container-interface.pbh>
#include <dune/xt/la/container/vector-interface.pbh>
#include <dune/xt/la/container/matrix-interface.pbh>
#include <dune/xt/la/solver.pbh>

#include <dune/xt/grid/walker/apply-on.bindings.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/gridprovider.pbh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/boundaryinfo.bindings.hh>
#include <dune/xt/grid/walker.bindings.hh>

#include <dune/xt/functions/interfaces.pbh>
#include <dune/xt/functions/expression.pbh>
#include <dune/xt/functions/constant.pbh>
#include <dune/xt/functions/checkerboard.pbh>
#include <dune/xt/functions/indicator.bindings.hh>
#include <dune/xt/functions/ESV2007.bindings.hh>
#include <dune/xt/functions/spe10.pbh>

#include <dune/gdt/operators/weighted-l2.bindings.hh>
#include <dune/gdt/operators/fluxreconstruction.bindings.hh>
#include <dune/gdt/operators/elliptic.bindings.hh>
#include <dune/gdt/operators/base.bindings.hh>
#include <dune/gdt/operators/elliptic-ipdg.bindings.hh>
#include <dune/gdt/operators/oswaldinterpolation.bindings.hh>
#include <dune/gdt/operators/l2.bindings.hh>
#include <dune/gdt/spaces/interface.bindings.hh>
#include <dune/gdt/spaces/constraints.bindings.hh>
#include <dune/gdt/prolongations.bindings.hh>
#include <dune/gdt/projections.bindings.hh>
#include <dune/gdt/projections/dirichlet.bindings.hh>
#include <dune/gdt/assembler/system.bindings.hh>
#include <dune/gdt/playground/spaces/block.bindings.hh>
#include <dune/gdt/functionals/base.bindings.hh>
#include <dune/gdt/functionals/elliptic-ipdg.bindings.hh>
#include <dune/gdt/functionals/l2.bindings.hh>
#include <dune/gdt/local/diffusive-flux-estimation-operator.bindings.hh>
#include <dune/gdt/local/elliptic-ipdg-operators.bindings.hh>
#include <dune/gdt/discretefunction/default.bindings.hh>


PYBIND11_PLUGIN(__usercode)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__usercode", "dune-gdt: user code");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.usercode");

  m.def("raise_SIGSEGV", []() { raise(SIGSEGV); });

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
