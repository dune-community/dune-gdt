// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/xt/grid/grids.hh>
#include <python/dune/xt/la/traits.hh>

#include "interfaces_all_grids.hh"


PYBIND11_MODULE(_operators_interfaces_istl_3d, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._spaces_interface");

//  OperatorInterface_for_all_grids<LA::IstlRowMajorSparseMatrix<double>,
//                                  LA::bindings::Istl,
//                                  void,
//                                  XT::Grid::bindings::Available3dGridTypes>::bind(m, "istl_sparse");
}
