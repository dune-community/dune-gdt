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
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/operators/l2.bindings.hh>
#include <dune/gdt/playground/spaces/restricted.hh>

using namespace Dune;
namespace py = pybind11;
using namespace pybind11::literals;
using Dune::XT::Grid::Layers;
using namespace Dune::XT;
using Dune::GDT::SpaceType;


template <class G,
          XT::Grid::Layers layer_type,
          XT::Grid::Backends layer_backend,
          size_t range_r = 1,
          size_t range_rC = 1,
          size_t source_r = range_r,
          size_t source_rC = range_rC>
void bind_l2_localizable_product(py::module& m)
{
  try {
    GDT::bindings::L2LocalizableProduct<G, layer_type, layer_backend, range_r, range_rC, source_r, source_rC>::bind(m);
  } catch (std::runtime_error&) {
  }
}


PYBIND11_PLUGIN(__operators_l2)
{
  py::module m("__operators_l2", "dune-gdt: L2LocalizableProduct, L2MatrixOperator");

  Dune::XT::Common::bindings::addbind_exceptions(m);

#if HAVE_DUNE_ALUGRID
  bind_l2_localizable_product<ALU_2D_SIMPLEX_CONFORMING, Layers::dd_subdomain, XT::Grid::Backends::part>(m);
  bind_l2_localizable_product<ALU_2D_SIMPLEX_CONFORMING, Layers::leaf, XT::Grid::Backends::part>(m);
#if HAVE_DUNE_FEM
  Dune::GDT::bindings::L2MatrixOperator<ALU_2D_SIMPLEX_CONFORMING,
                                        Layers::dd_subdomain,
                                        SpaceType::block_dg,
                                        GDT::Backends::fem,
                                        1,
                                        1,
                                        LA::Backends::istl_sparse>::bind(m);
  Dune::GDT::bindings::L2MatrixOperator<ALU_2D_SIMPLEX_CONFORMING,
                                        Layers::dd_subdomain,
                                        SpaceType::dg,
                                        GDT::Backends::fem,
                                        1,
                                        1,
                                        LA::Backends::istl_sparse>::bind(m);
  Dune::GDT::bindings::L2MatrixOperator<ALU_2D_SIMPLEX_CONFORMING,
                                        Layers::leaf,
                                        SpaceType::dg,
                                        GDT::Backends::fem,
                                        1,
                                        1,
                                        LA::Backends::istl_sparse>::bind(m);
  Dune::GDT::bindings::L2MatrixOperator<ALU_2D_SIMPLEX_CONFORMING,
                                        Layers::level,
                                        SpaceType::dg,
                                        GDT::Backends::fem,
                                        1,
                                        1,
                                        LA::Backends::istl_sparse>::bind(m);
  Dune::GDT::bindings::internal::
      L2MatrixOperator<GDT::RestrictedSpace<
                           typename GDT::SpaceProvider<ALU_2D_SIMPLEX_CONFORMING,
                                                       Layers::leaf,
                                                       GDT::SpaceType::rt,
                                                       GDT::Backends::pdelab,
                                                       0,
                                                       double,
                                                       2>::type,
                           typename XT::Grid::Layer<ALU_2D_SIMPLEX_CONFORMING,
                                                    Layers::dd_subdomain,
                                                    XT::Grid::Backends::part,
                                                    XT::Grid::DD::SubdomainGrid<ALU_2D_SIMPLEX_CONFORMING>>::type>,
                       XT::LA::IstlRowMajorSparseMatrix<double>>::
          bind(m, "RtPdelabAlu2dSimplexLeafRestrictedSubdomainPartSpace", "istl_row_major_sparse_matrix_double");
#endif // HAVE_DUNE_FEM
#endif // HAVE_DUNE_ALUGRID

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = Dune::XT::Common::numeric_cast<int>(args.size());
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
        "max_info_level"_a = std::numeric_limits<ssize_t>::max(),
        "max_debug_level"_a = std::numeric_limits<ssize_t>::max(),
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  m.def("_test_logger",
        [](const bool info, const bool debug, const bool warning) {
          auto logger = Dune::XT::Common::TimedLogger().get("dune.gdt.operators.elliptic");
          if (info)
            logger.info() << "info logging works!" << std::endl;
          if (debug)
            logger.debug() << "debug logging works!" << std::endl;
          if (warning)
            logger.warn() << "warning logging works!" << std::endl;
        },
        "info"_a = true,
        "debug"_a = true,
        "warning"_a = true);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
