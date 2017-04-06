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

#include <dune/xt/common/exceptions.hh>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/gdt/operators/elliptic-ipdg.bindings.hh>


PYBIND11_PLUGIN(__operators_elliptic_ipdg)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__operators_elliptic_ipdg", "dune-gdt: EllipticIpdgMatrixOperator");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

// alu_fem_istl.cc
#if HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, leaf, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, level, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, dd_subdomain, part, cg, fem, 1, istl_sparse);
#endif

// yasp_fem_istl.cc
#if HAVE_DUNE_FEM && HAVE_DUNE_ISTL
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_YASP(m, leaf, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_YASP(m, level, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_YASP(m, dd_subdomain, part, cg, fem, 1, istl_sparse);
#endif

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = boost::numeric_cast<int>(args.size());
          char** argv = Dune::XT::Common::vector_to_main_args(args);
          Dune::MPIHelper::instance(argc, argv);
#if HAVE_DUNE_FEM
          Dune::Fem::MPIManager::initialize(argc, argv);
#endif
        },
        "args"_a = std::vector<std::string>());

  m.def("_init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = -1,
        "max_debug_level"_a = -1,
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
