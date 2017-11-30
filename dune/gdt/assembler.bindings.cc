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

#include <dune/gdt/assembler/system.bindings.hh>
#include <dune/gdt/spaces/constraints.bindings.hh>


#define DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, _GRID, _layer, _backend, _layer_name)                                     \
  auto dirichlet_constraints_##_GRID##_##_layer##_##_backend = Dune::GDT::bindings::                                   \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::_backend,                      \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::bind(_m, _layer_name)

#define DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_GRID, _layer, _backend, _la)                                           \
  Dune::GDT::bindings::                                                                                                \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::_backend,                      \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::                                                                                    \
          addbind<Dune::XT::LA::Backends::_la>(dirichlet_constraints_##_GRID##_##_layer##_##_backend)


PYBIND11_PLUGIN(__assembler)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__assembler", "dune-gdt: SystemAssembler");

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::class_<Dune::GDT::bindings::ResultStorage> ResultStorage(m, "ResultStorage", "dune-gdt: ResultStorage");
  ResultStorage.def(pybind11::init<>());
  ResultStorage.def_property(
      "result",
      [](const Dune::GDT::bindings::ResultStorage& self) { return self.result(); },
      [](Dune::GDT::bindings::ResultStorage& self, const double& value) { self.result() = value; });

  DUNE_GDT_SPACES_CONSTRAINTS_BIND(m, ALU_2D_SIMPLEX_CONFORMING, leaf, view, "leaf");
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(ALU_2D_SIMPLEX_CONFORMING, leaf, view, istl_sparse);
  DUNE_GDT_SPACES_CONSTRAINTS_BIND(m, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, part, "dd_subdomain");
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, part, istl_sparse);

  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(m);

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
          auto logger = Dune::XT::Common::TimedLogger().get("dune.gdt.assembler");
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
