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
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "binary_intersection_combined.hh"
#include "binary_intersection_interface.hh"


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct LocalBinaryIntersectionIntegrandInterface_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;
  using F = double;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::LocalBinaryIntersectionIntegrandInterface;
    using Dune::GDT::bindings::LocalBinaryIntersectionIntegrandSum;
    using Dune::GDT::bindings::LocalBinaryToUnaryIntersectionIntegrand;

    LocalBinaryIntersectionIntegrandInterface<G, I>::bind(m, "leaf");
    if (d > 1) {
      LocalBinaryIntersectionIntegrandInterface<G, I, 1, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, 1, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandInterface<G, I, d, d, F, F, d, d, F>::bind(m, "leaf");
    }
    // add your extra dimensions here
    // ...
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CI = typename GridGlueType::GlueType::Intersection;
      using CCI = Dune::XT::Grid::internal::CouplingIntersectionWithCorrectNormal<CI, I>;
      LocalBinaryIntersectionIntegrandInterface<G, CCI>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, 1, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, 1, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, 1, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, d, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, d, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandInterface<G, CCI, d, d, F, F, d, d, F>::bind(m, "coupling");
    }
#endif

    // we need to bind LocalBinaryToUnaryIntersectionIntegrand here, since it requires the above interface
    LocalBinaryToUnaryIntersectionIntegrand<G, I>::bind(m, "leaf");
    if (d > 1) {
      LocalBinaryToUnaryIntersectionIntegrand<G, I, 1, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, 1, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryToUnaryIntersectionIntegrand<G, I, d, d, F, F, d, d, F>::bind(m, "leaf");
    }
    // add your extra dimensions here
    // ...
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CI = typename GridGlueType::GlueType::Intersection;
      using CCI = Dune::XT::Grid::internal::CouplingIntersectionWithCorrectNormal<CI, I>;
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, 1, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, 1, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, 1, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, d, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, d, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryToUnaryIntersectionIntegrand<G, CCI, d, d, F, F, d, d, F>::bind(m, "coupling");
    }
#endif
    // we need to bind LocalBinaryIntersectionIntegrandSum here, since it requires the above interface
    LocalBinaryIntersectionIntegrandSum<G, I>::bind(m, "leaf");
    if (d > 1) {
      LocalBinaryIntersectionIntegrandSum<G, I, 1, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, 1, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, 1, F, F, d, d, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, 1, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, d, 1, F>::bind(m, "leaf");
      LocalBinaryIntersectionIntegrandSum<G, I, d, d, F, F, d, d, F>::bind(m, "leaf");
    }
    // add your extra dimensions here
    // ...
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CI = typename GridGlueType::GlueType::Intersection;
      using CCI = Dune::XT::Grid::internal::CouplingIntersectionWithCorrectNormal<CI, I>;
      LocalBinaryIntersectionIntegrandSum<G, CCI>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, 1, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, 1, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, 1, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, 1, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, 1, F, F, d, d, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, d, F, F, 1, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, d, F, F, d, 1, F>::bind(m, "coupling");
      LocalBinaryIntersectionIntegrandSum<G, CCI, d, d, F, F, d, d, F>::bind(m, "coupling");
    }
#endif
    LocalBinaryIntersectionIntegrandInterface_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct LocalBinaryIntersectionIntegrandInterface_for_all_grids<Dune::XT::Common::tuple_null_type>
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

  LocalBinaryIntersectionIntegrandInterface_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
