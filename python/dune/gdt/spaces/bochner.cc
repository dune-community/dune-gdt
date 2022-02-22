// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/spaces/bochner.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


/**
 * \note Assumes that GV is the leaf view!
 */
template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class BochnerSpace
{
  using G = typename GV::Grid;
  static const size_t d = G::dimension;
  using SpatialSpaceType = GDT::SpaceInterface<GV, r, rC, R>;

public:
  using type = GDT::BochnerSpace<GV, r, rC, R>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "bochner_space")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    if (r > 1)
      class_name += "_to_" + XT::Common::to_string(r) + "d";
    if (!std::is_same<R, double>::value)
      class_name += "_" + XT::Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](const SpatialSpaceType& spatial_space, const std::vector<double>& time_points) {
            return new type(spatial_space, time_points);
          }),
          "spatial_space"_a,
          "time_points"_a,
          py::keep_alive<1, 2>()); // spatial_space
    c.def_property_readonly("spatial_space", &type::spatial_space);
    c.def_property_readonly("temporal_space", &type::temporal_space);
    c.def_property_readonly("time_interval", &type::time_interval);

    std::string space_type_name = class_id;
    if (!std::is_same<R, double>::value)
      space_type_name += "_" + XT::Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    m.def(
        XT::Common::to_camel_case(space_type_name).c_str(),
        [](const SpatialSpaceType& spatial_space, const std::vector<double>& time_points) {
          return new type(spatial_space, time_points);
        },
        "spatial_space"_a,
        "time_points"_a,
        py::keep_alive<0, 1>()); // spatial_space

    return c;
  } // ... bind(...)
}; // class BochnerSpace


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct BochnerSpace_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::BochnerSpace;
    using Dune::XT::Grid::bindings::grid_name;

    BochnerSpace<GV>::bind(m, grid_name<G>::value());
    // add your extra dimensions here
    // ...
    BochnerSpace_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct BochnerSpace_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_spaces_bochner, m)
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

  BochnerSpace_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
