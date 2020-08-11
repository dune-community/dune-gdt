// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
#define PYTHON_DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

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
template <class GV, size_t r = 1, class R = double>
class ContinuousLagrangeSpace
{
  using G = typename GV::Grid;
  static const size_t d = G::dimension;

public:
  using type = GDT::ContinuousLagrangeSpace<GV, r, R>;
  using base_type = GDT::SpaceInterface<GV, r, 1, R>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type
  bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "continuous_lagrange_space")
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
    c.def(py::init([](XT::Grid::GridProvider<G>& grid_provider, const int order, const std::string& logging_prefix) {
            return new type(grid_provider.leaf_view(), order, logging_prefix);
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
    if (r == 1)
      m.def(
          XT::Common::to_camel_case(space_type_name).c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const int order,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.leaf_view(), order, logging_prefix); },
          "grid"_a,
          "order"_a,
          "dim_range"_a = XT::Grid::bindings::Dimension<r>(),
          "logging_prefix"_a = "");
    else
      m.def(
          XT::Common::to_camel_case(space_type_name).c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const int order,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.leaf_view(), order, logging_prefix); },
          "grid"_a,
          "order"_a,
          "dim_range"_a,
          "logging_prefix"_a = "");

    return c;
  } // ... bind(...)
}; // class ContinuousLagrangeSpace


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_SPACES_H1_CONTINUOUS_LAGRANGE_HH
