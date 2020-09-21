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

#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/spaces/hdiv/raviart-thomas.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>

namespace Dune {
namespace GDT {
namespace bindings {


/**
 * \note Assumes that GV is the leaf view!
 */
template <class GV, class R = double>
class RaviartThomasSpace
{
  using G = typename GV::Grid;
  static const size_t d = G::dimension;

public:
  using type = GDT::RaviartThomasSpace<GV, R>;
  using base_type = GDT::SpaceInterface<GV, d, 1, R>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type
  bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "raviart_thomas_space")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    class_name += "_to_" + XT::Common::to_string(size_t(d)) + "d";
    if (!std::is_same<R, double>::value)
      class_name += "_" + XT::Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](XT::Grid::GridProvider<G>& grid_provider, const int order, const std::string& logging_prefix) {
            return new type(grid_provider.leaf_view(), order, logging_prefix); // Otherwise we get an error here!
          }),
          "grid_provider"_a,
          "order"_a,
          "logging_prefix"_a = "");
    c.def("__repr__", [](const type& self) {
      std::stringstream ss;
      ss << self;
      return ss.str();
    });

    std::string space_type_name = class_id;
    if (!std::is_same<R, double>::value)
      space_type_name += "_" + XT::Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    m.def(
        XT::Common::to_camel_case(space_type_name).c_str(),
        [](XT::Grid::GridProvider<G>& grid, const int order, const std::string& logging_prefix) {
          return new type(grid.leaf_view(), order, logging_prefix); // Otherwise we get an error here!
        },
        "grid"_a,
        "order"_a,
        "logging_prefix"_a = "");

    return c;
  } // ... bind(...)
}; // class RaviartThomasSpace


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct RaviartThomasSpace_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::RaviartThomasSpace;
    using Dune::XT::Grid::bindings::grid_name;

    RaviartThomasSpace<GV>::bind(m, grid_name<G>::value());

    RaviartThomasSpace_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct RaviartThomasSpace_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_spaces_hdiv_raviart_thomas, m)
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

  RaviartThomasSpace_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
