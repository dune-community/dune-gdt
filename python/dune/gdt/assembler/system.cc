// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/available_types.hh>
#include <python/dune/gdt/shared.hh>

#include <python/dune/gdt/assembler/system.hh>
#include <python/dune/gdt/spaces/constraints.hh>


template <class _GRID, Dune::XT::Grid::Layers _layer, Dune::XT::Grid::Backends _backend, Dune::XT::LA::Backends _la>
void bind_constraints(pybind11::module& _m, std::string _layer_name)
{
  using Intersection = Dune::XT::Grid::extract_intersection_t<
      typename Dune::XT::Grid::Layer<_GRID, _layer, _backend, Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>;
  using DC = Dune::GDT::bindings::DirichletConstraints<Intersection, _GRID>;
  auto dc = DC::bind(_m, _layer_name);
  DC::template addbind<_la>(dc);
}

template <class _G,
          Dune::XT::Grid::Layers _gl,
          Dune::XT::Grid::Backends _gl_backend,
          Dune::GDT::Backends _t_backend,
          Dune::GDT::SpaceType _t_type,
          Dune::XT::Grid::Layers _t_gl,
          size_t _t_p,
          size_t _t_r,
          Dune::GDT::Backends _a_backend,
          Dune::GDT::SpaceType _a_type,
          Dune::XT::Grid::Layers _a_gl,
          size_t _a_p,
          size_t _a_r>
void bind_system(pybind11::module& _m)
{
  using TestSpc = Dune::GDT::SpaceProvider<_G, _t_gl, _t_type, _t_backend, _t_p, double, _t_r, 1>;
  using AnsatzSpc = Dune::GDT::SpaceProvider<_G, _a_gl, _a_type, _a_backend, _a_p, double, _a_r, 1>;
  Dune::GDT::bindings::SystemAssembler<TestSpc, _gl, _gl_backend, AnsatzSpc>::bind(_m);
}

template <class Tuple>
void addbind_for_Grid(pybind11::module& m)
{
  using namespace Dune::XT::Grid;
  using G = typename Tuple::head_type;
  constexpr auto leaf = Dune::XT::Grid::Layers::leaf;
  constexpr auto dd_subdomain = Dune::XT::Grid::Layers::dd_subdomain;
  constexpr auto dd_subdomain_oversampled = Dune::XT::Grid::Layers::dd_subdomain_oversampled;
  constexpr auto dd_subdomain_boundary = Dune::XT::Grid::Layers::dd_subdomain_boundary;
  constexpr auto dd_subdomain_coupling = Dune::XT::Grid::Layers::dd_subdomain_coupling;
  constexpr auto view = Dune::XT::Grid::Backends::view;
  constexpr auto istl_sparse = Dune::XT::LA::Backends::istl_sparse;
  constexpr auto gdt = Dune::GDT::Backends::gdt;
  constexpr auto dg = Dune::GDT::SpaceType::dg;

  bind_constraints<G, leaf, view, istl_sparse>(m, "leaf");
  bind_constraints<G, dd_subdomain, view, istl_sparse>(m, "dd_subdomain");

  bind_system<G, leaf, view, gdt, dg, leaf, 1, 1, gdt, dg, leaf, 1, 1>(m);
  bind_system<G, leaf, view, gdt, dg, leaf, 1, 1, gdt, dg, leaf, 1, 1>(m);
  bind_system<G, leaf, view, gdt, dg, leaf, 2, 1, gdt, dg, leaf, 2, 1>(m);
  bind_system<G, leaf, view, gdt, dg, leaf, 3, 1, gdt, dg, leaf, 3, 1>(m);
  bind_system<G, dd_subdomain, view, gdt, dg, dd_subdomain, 1, 1, gdt, dg, dd_subdomain, 1, 1>(m);
  bind_system<G, dd_subdomain_boundary, view, gdt, dg, dd_subdomain, 1, 1, gdt, dg, dd_subdomain, 1, 1>(m);
  bind_system<G, dd_subdomain_coupling, view, gdt, dg, dd_subdomain, 1, 1, gdt, dg, dd_subdomain, 1, 1>(m);
  bind_system<G, dd_subdomain_oversampled, view, gdt, dg, dd_subdomain, 1, 1, gdt, dg, dd_subdomain, 1, 1>(m);

  addbind_for_Grid<typename Tuple::tail_type>(m);
} // ... addbind_for_Grid(...)

template <>
void addbind_for_Grid<boost::tuples::null_type>(pybind11::module&)
{}

PYBIND11_MODULE(__assembler, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.assembler");

  py::class_<Dune::GDT::bindings::ResultStorage> ResultStorage(m, "ResultStorage", "dune-gdt: ResultStorage");
  ResultStorage.def(pybind11::init<>());
  ResultStorage.def_property(
      "result",
      [](const Dune::GDT::bindings::ResultStorage& self) { return self.result(); },
      [](Dune::GDT::bindings::ResultStorage& self, const double& value) { self.result() = value; });

  py::module::import("dune.xt.grid.walker");
  using TP = boost::tuple<GDT_BINDINGS_GRID>;
  addbind_for_Grid<TP>(m);
}
