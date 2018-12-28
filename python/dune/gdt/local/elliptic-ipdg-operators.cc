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
#include <python/dune/gdt/shared.hh>

#include "elliptic-ipdg-operators.hh"

template <class _G,
          size_t _d,
          Dune::XT::Grid::Layers _g_layer,
          Dune::XT::Grid::Backends _g_backend,
          Dune::GDT::SpaceType _s_type,
          Dune::GDT::Backends _s_backend,
          Dune::XT::Grid::Layers _i_layer_type,
          int _p,
          class _R,
          size_t _r,
          size_t _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method _method>
void bind_1d(pybind11::module& _m, std::string _layer_name)
{
  using F1 = Dune::XT::Functions::LocalizableFunctionInterface<
      Dune::XT::Grid::extract_entity_t<
          typename Dune::XT::Grid::Layer<_G, _g_layer, _g_backend, Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>,
      double,
      _d,
      double,
      1,
      1>;
  using Fd = Dune::XT::Functions::LocalizableFunctionInterface<
      Dune::XT::Grid::extract_entity_t<
          typename Dune::XT::Grid::Layer<_G, _g_layer, _g_backend, Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>,
      double,
      _d,
      double,
      _d,
      _d>;
  Dune::GDT::bindings::LocalEllipticIpdgInnerIntegralOperator<
      F1,
      Fd,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
  Dune::GDT::bindings::LocalEllipticIpdgInnerIntegralOperator<
      F1,
      void,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
  Dune::GDT::bindings::LocalEllipticIpdgBoundaryIntegralOperator<
      F1,
      Fd,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
  Dune::GDT::bindings::LocalEllipticIpdgBoundaryIntegralOperator<
      F1,
      void,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
}

template <class _G,
          size_t _d,
          Dune::XT::Grid::Layers _g_layer,
          Dune::XT::Grid::Backends _g_backend,
          Dune::GDT::SpaceType _s_type,
          Dune::GDT::Backends _s_backend,
          Dune::XT::Grid::Layers _i_layer_type,
          int _p,
          class _R,
          size_t _r,
          size_t _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method _method>
void bind_dd(pybind11::module& _m, std::string _layer_name)
{
  bind_1d<_G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _method>(_m, _layer_name);
  using Function = Dune::XT::Functions::LocalizableFunctionInterface<
      Dune::XT::Grid::extract_entity_t<
          typename Dune::XT::Grid::Layer<_G, _g_layer, _g_backend, Dune::XT::Grid::DD::SubdomainGrid<_G>>::type>,
      double,
      _d,
      double,
      _d,
      _d>;
  Dune::GDT::bindings::LocalEllipticIpdgInnerIntegralOperator<
      Function,
      void,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
  Dune::GDT::bindings::LocalEllipticIpdgBoundaryIntegralOperator<
      Function,
      void,
      Dune::GDT::SpaceProvider<_G, _g_layer, _s_type, _s_backend, _p, _R, _r, _rC>,
      _i_layer_type,
      _method>::bind(_m, _layer_name);
}


template <class _G,
          size_t _d,
          Dune::XT::Grid::Layers _g_layer,
          Dune::XT::Grid::Backends _g_backend,
          Dune::GDT::SpaceType _s_type,
          Dune::GDT::Backends _s_backend,
          Dune::XT::Grid::Layers _i_layer_type,
          int _p,
          class _R,
          size_t _r,
          size_t _rC>
void bind_ipdg(pybind11::module& _m, std::string _layer_name)
{
  bind_dd<_G,
          _d,
          _g_layer,
          _g_backend,
          _s_type,
          _s_backend,
          _i_layer_type,
          _p,
          _R,
          _r,
          _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method::sipdg>(_m, _layer_name);
  bind_dd<_G,
          _d,
          _g_layer,
          _g_backend,
          _s_type,
          _s_backend,
          _i_layer_type,
          _p,
          _R,
          _r,
          _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method::swipdg>(_m, _layer_name);
  bind_dd<_G,
          _d,
          _g_layer,
          _g_backend,
          _s_type,
          _s_backend,
          _i_layer_type,
          _p,
          _R,
          _r,
          _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(_m, _layer_name);
  bind_dd<_G,
          _d,
          _g_layer,
          _g_backend,
          _s_type,
          _s_backend,
          _i_layer_type,
          _p,
          _R,
          _r,
          _rC,
          Dune::GDT::LocalEllipticIpdgIntegrands::Method::swipdg_affine_tensor>(_m, _layer_name);
}


PYBIND11_MODULE(__local_elliptic_ipdg_operators, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.grid.provider");
  py::module::import("dune.gdt.__spaces");
  py::module::import("dune.gdt.__spaces_block");
  bind_ipdg<GDT_BINDINGS_GRID,
            2,
            Dune::XT::Grid::Layers::dd_subdomain,
            Dune::XT::Grid::Backends::view,
            Dune::GDT::SpaceType::dg,
            Dune::GDT::Backends::gdt,
            Dune::XT::Grid::Layers::dd_subdomain_coupling,
            1,
            double,
            1,
            1>(m, "_dd_subdomain_");

  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.assembler");
}
