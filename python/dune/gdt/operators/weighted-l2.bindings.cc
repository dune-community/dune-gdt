// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/operators/weighted-l2.bindings.hh>


PYBIND11_MODULE(__operators_weighted_l2, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using Dune::XT::Grid::Layers;
  using Dune::XT::Grid::Backends;


#if HAVE_DUNE_ALUGRID
  Dune::GDT::bindings::WeightedL2LocalizableProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, Backends::view>::bind(m);
  Dune::GDT::bindings::WeightedL2LocalizableProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::level, Backends::view>::bind(m);
  Dune::GDT::bindings::WeightedL2LocalizableProduct<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, Backends::view>::
      bind(m);
#endif // HAVE_DUNE_ALUGRID
}

#endif // HAVE_DUNE_PYBINDXI
