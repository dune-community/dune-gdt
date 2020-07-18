// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "binary_intersection_combined.hh"
#include "binary_intersection_interface.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalBinaryIntersectionIntegrandInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;
  using F = double;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::LocalBinaryIntersectionIntegrandInterface;
    using Dune::GDT::bindings::LocalBinaryIntersectionIntegrandSum;
    using Dune::GDT::bindings::LocalBinaryToUnaryIntersectionIntegrand;

    LocalBinaryIntersectionIntegrandInterface<G, I>::bind(m);
    if (d > 1) {
      LocalBinaryIntersectionIntegrandInterface<G, I, 1, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, 1, 1, F, F, d, d, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, 1, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, d, d, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, 1, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, d, d, F>::bind(m);
    }
    // add your extra dimensions here
    // ...

    // we need to bind LocalBinaryToUnaryIntersectionIntegrand here, since it requires the above interface
    LocalBinaryToUnaryIntersectionIntegrand<G, I>::bind(m);
    if (d > 1) {
      LocalBinaryToUnaryIntersectionIntegrand<G, I, 1, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, 1, 1, F, F, d, d, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, 1, 1, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, d, d, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, 1, 1, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, d, 1, F>::bind(m);
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, d, d, F>::bind(m);
    }
    // add your extra dimensions here
    // ...

    // we need to bind LocalBinaryIntersectionIntegrandSum here, since it requires the above interface
    LocalBinaryIntersectionIntegrandSum<G, I>::bind(m);
    if (d > 1) {
      LocalBinaryIntersectionIntegrandSum<G, I, 1, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, 1, 1, F, F, d, d, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, 1, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, d, d, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, 1, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, d, 1, F>::bind(m);
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, d, d, F>::bind(m);
    }
    // add your extra dimensions here
    // ...

    LocalBinaryIntersectionIntegrandInterface_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalBinaryIntersectionIntegrandInterface_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_binary_intersection_interface, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_integrands_unary_intersection_interface");

  LocalBinaryIntersectionIntegrandInterface_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
