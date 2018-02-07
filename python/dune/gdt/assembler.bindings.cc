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
#include <python/dune/gdt/shared.hh>

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


#define DUNE_GDT_ASSEMBLER_SYSTEM_BIND(                                                                                \
    _m, _G, _gl, _gl_backend, _t_backend, _t_type, _t_gl, _t_p, _t_r, _a_backend, _a_type, _a_gl, _a_p, _a_r)          \
  Dune::GDT::bindings::SystemAssembler<Dune::GDT::SpaceProvider<_G,                                                    \
                                                                Dune::XT::Grid::Layers::_t_gl,                         \
                                                                Dune::GDT::SpaceType::_t_type,                         \
                                                                Dune::GDT::Backends::_t_backend,                       \
                                                                _t_p,                                                  \
                                                                double,                                                \
                                                                _t_r,                                                  \
                                                                1>,                                                    \
                                       Dune::XT::Grid::Layers::_gl,                                                    \
                                       Dune::XT::Grid::Backends::_gl_backend,                                          \
                                       Dune::GDT::SpaceProvider<_G,                                                    \
                                                                Dune::XT::Grid::Layers::_a_gl,                         \
                                                                Dune::GDT::SpaceType::_a_type,                         \
                                                                Dune::GDT::Backends::_a_backend,                       \
                                                                _a_p,                                                  \
                                                                double,                                                \
                                                                _a_r,                                                  \
                                                                1>>::bind(_m)

PYBIND11_MODULE(__assembler, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::class_<Dune::GDT::bindings::ResultStorage> ResultStorage(m, "ResultStorage", "dune-gdt: ResultStorage");
  ResultStorage.def(pybind11::init<>());
  ResultStorage.def_property(
      "result",
      [](const Dune::GDT::bindings::ResultStorage& self) { return self.result(); },
      [](Dune::GDT::bindings::ResultStorage& self, const double& value) { self.result() = value; });

  using G = ALU_2D_SIMPLEX_CONFORMING;
  add_initialization(m, "dune.gdt.assembler");
  DUNE_GDT_SPACES_CONSTRAINTS_BIND(m, G, leaf, view, "leaf");
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(G, leaf, view, istl_sparse);
  DUNE_GDT_SPACES_CONSTRAINTS_BIND(m, G, dd_subdomain, part, "dd_subdomain");
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(G, dd_subdomain, part, istl_sparse);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(m, G, leaf, part, fem, dg, leaf, 1, 1, fem, dg, leaf, 1, 1);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(m, G, leaf, part, fem, dg, leaf, 2, 1, fem, dg, leaf, 2, 1);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(m, G, leaf, part, fem, dg, leaf, 3, 1, fem, dg, leaf, 3, 1);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(m, G, dd_subdomain, part, fem, dg, dd_subdomain, 1, 1, fem, dg, dd_subdomain, 1, 1);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(
      m, G, dd_subdomain_boundary, part, fem, dg, dd_subdomain, 1, 1, fem, dg, dd_subdomain, 1, 1);
  DUNE_GDT_ASSEMBLER_SYSTEM_BIND(
      m, G, dd_subdomain_coupling, part, fem, dg, dd_subdomain, 1, 1, fem, dg, dd_subdomain, 1, 1);
}

#endif // HAVE_DUNE_PYBINDXI
