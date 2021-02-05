// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <python/dune/xt/grid/grids.bindings.hh>

#include "matrix-based_factory.hh"


PYBIND11_MODULE(_operators_matrix_based_factory, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_bilinear_forms_element_interface");
  py::module::import("dune.gdt._operators_interfaces_common");
  py::module::import("dune.gdt._operators_interfaces_eigen");
  py::module::import("dune.gdt._operators_interfaces_istl_1d");
  py::module::import("dune.gdt._operators_interfaces_istl_2d");
  py::module::import("dune.gdt._operators_interfaces_istl_3d");
  py::module::import("dune.gdt._spaces_interface");

  //  MatrixOperatorFactory_for_all_grids<LA::CommonDenseMatrix<double>,
  //                               LA::bindings::Common,
  //                               LA::bindings::Dense,
  //                               XT::Grid::bindings::AvailableGridTypes>::bind(m, "common_dense");
  //  // Generic linear solver missing for CommonSparseMatrix!
  //  MatrixOperatorFactory_for_all_grids<LA::CommonSparseMatrix<double>,
  //                               LA::bindings::Common,
  //                               LA::bindings::Sparse,
  //                               XT::Grid::bindings::AvailableGridTypes>::bind(m, "common_sparse");
  //#if HAVE_EIGEN
  //  MatrixOperatorFactory_for_all_grids<LA::EigenDenseMatrix<double>,
  //                               LA::bindings::Eigen,
  //                               LA::bindings::Dense,
  //                               XT::Grid::bindings::AvailableGridTypes>::bind(m, "eigen_dense");
  //  MatrixOperatorFactory_for_all_grids<LA::EigenRowMajorSparseMatrix<double>,
  //                               LA::bindings::Eigen,
  //                               LA::bindings::Sparse,
  //                               XT::Grid::bindings::AvailableGridTypes>::bind(m, "eigen_sparse");
  //#endif
//  MatrixOperatorFactory_for_all_grids<LA::IstlRowMajorSparseMatrix<double>,
//                                      LA::bindings::Istl,
//                                      void,
//                                      XT::Grid::bindings::AvailableGridTypes>::bind(m, "istl_sparse");
}
