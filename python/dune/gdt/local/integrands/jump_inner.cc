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
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/local/integrands/jump.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>

#include "quaternary_intersection_interface.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class G, class I, size_t r = 1>
class LocalInnerJumpIntegrand
{
  using E = XT::Grid::extract_entity_t<G>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalJumpIntegrands::Inner<I, r>;
  using base_type = typename type::BaseType;
  using bound_type = pybind11::class_<type, base_type>;
  using F = typename type::F;

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "local_inner_jump_integrand")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto class_name =
        class_id + "_" + LocalQuaternaryIntersectionIntegrandInterface<G, I, r>::id(grid_id, layer_id);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](const bool weighted_by_intersection_diameter, const std::string& logging_prefix) {
            if (weighted_by_intersection_diameter)
              return new type(GDT::LocalIPDGIntegrands::internal::default_intersection_diameter<I>(), logging_prefix);
            else
              return new type([](const auto&) { return 1.; }, logging_prefix);
          }),
          "weighted_by_intersection_diameter"_a = true,
          "logging_prefix"_a = "");

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    if constexpr (r == 1)
      m.def(
          FactoryName.c_str(),
          [](const GP& /*grid*/,
             const XT::Grid::bindings::Dimension<r>&,
             const bool weighted_by_intersection_diameter,
             const std::string& logging_prefix) {
            if (weighted_by_intersection_diameter)
              return new type(GDT::LocalIPDGIntegrands::internal::default_intersection_diameter<I>(), logging_prefix);
            else
              return new type([](const auto&) { return 1.; }, logging_prefix);
          },
          "grid"_a,
          "dim_range"_a = XT::Grid::bindings::Dimension<r>(),
          "weighted_by_intersection_diameter"_a = true,
          "logging_prefix"_a = "");
    else
      m.def(
          FactoryName.c_str(),
          [](const GP&,
             const XT::Grid::bindings::Dimension<r>&,
             const bool weighted_by_intersection_diameter,
             const std::string& logging_prefix) {
            if (weighted_by_intersection_diameter)
              return new type(GDT::LocalIPDGIntegrands::internal::default_intersection_diameter<I>(), logging_prefix);
            else
              return new type([](const auto&) { return 1.; }, logging_prefix);
          },
          "grid"_a,
          "dim_range"_a,
          "weighted_by_intersection_diameter"_a = true,
          "logging_prefix"_a = "");

    return c;
  } // ... bind(...)
}; // class LocalInnerJumpIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct LocalInnerJumpIntegrand_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalInnerJumpIntegrand<G, I>::bind(m);

    if (d > 1)
      Dune::GDT::bindings::LocalInnerJumpIntegrand<G, I, d>::bind(m);
    // add your extra dimensions here

    LocalInnerJumpIntegrand_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct LocalInnerJumpIntegrand_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_jump_inner, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt._local_integrands_quaternary_intersection_interface");

  LocalInnerJumpIntegrand_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
