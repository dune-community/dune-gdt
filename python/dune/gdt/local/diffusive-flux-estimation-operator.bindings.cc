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

#include "diffusive-flux-estimation-operator.bindings.hh"


#define DUNE_GDT_LOCAL_DIFFUSIVE_FLUX_ESTIMATION_OPERATOR_BIND(_G, _m)                                                 \
  Dune::GDT::bindings::LocalDiffusiveFluxEstimationOperator<_G>::bind(_m)

PYBIND11_MODULE(__local_diffusive_flux_estimation_operator)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.local.operators.diffusive-flux-reconstruction");
  DUNE_GDT_LOCAL_DIFFUSIVE_FLUX_ESTIMATION_OPERATOR_BIND(ALU_2D_SIMPLEX_CONFORMING, m);
}

#endif // HAVE_DUNE_PYBINDXI
