// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH
#define PYTHON_DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

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
class DiscontinuousLagrangeSpace
{
  using G = typename GV::Grid;
  static const size_t d = G::dimension;

public:
  using type = GDT::DiscontinuousLagrangeSpace<GV, r, R>;
  using base_type = GDT::SpaceInterface<GV, r, 1, R>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type
  bind(pybind11::module& m, const std::string& grid_id, const std::string& class_id = "discontinuous_lagrange_space")
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
    c.def_property_readonly("dimwise_global_mapping", [](type& self) { return self.dimwise_global_mapping; });
    c.def(py::init([](XT::Grid::GridProvider<G>& grid_provider, const int order, const bool dimwise_global_mapping) {
            return new type(grid_provider.leaf_view(), order, dimwise_global_mapping); // Otherwise we get an error
                                                                                       // here!
          }),
          "grid_provider"_a,
          "order"_a,
          "dimwise_global_mapping"_a = (r == 1) ? false : true);
    c.def(
        "interpolation_points",
        [](type& self) {
          DUNE_THROW_IF(self.min_polorder() != self.max_polorder(),
                        XT::Common::Exceptions::wrong_input_given,
                        "Not implemented for varying polynomial degree!");
          const GDT::DiscontinuousLagrangeSpace<GV, 1, R> scalar_self(self.grid_view(), self.min_polorder());
          const auto& global_basis = scalar_self.basis();
          auto basis = global_basis.localize();
          DynamicVector<size_t> global_DoF_indices(scalar_self.mapper().max_local_size());
          auto points = std::make_unique<XT::LA::CommonDenseMatrix<double>>(scalar_self.mapper().size(), size_t(d), 0.);
          for (auto&& element : elements(scalar_self.grid_view())) {
            basis->bind(element);
            scalar_self.mapper().global_indices(element, global_DoF_indices);
            auto local_lagrange_points = basis->finite_element().lagrange_points();
            for (size_t ii = 0; ii < basis->size(); ++ii) {
              auto global_point = element.geometry().global(local_lagrange_points[ii]);
              for (size_t jj = 0; jj < size_t(d); ++jj)
                points->set_entry(global_DoF_indices[ii], jj, global_point[jj]);
            }
          }
          return std::move(points);
        },
        py::call_guard<py::gil_scoped_release>());
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
             const bool dimwise_global_mapping,
             const XT::Grid::bindings::Dimension<r>&) {
            return new type(grid.leaf_view(), order, dimwise_global_mapping); // Otherwise we get an error here!
          },
          "grid"_a,
          "order"_a,
          "dimwise_global_mapping"_a = (r == 1) ? false : true,
          "dim_range"_a = XT::Grid::bindings::Dimension<r>());
    else
      m.def(
          XT::Common::to_camel_case(space_type_name).c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const int order,
             const bool dimwise_global_mapping,
             const XT::Grid::bindings::Dimension<r>&) {
            return new type(grid.leaf_view(), order, dimwise_global_mapping); // Otherwise we get an error here!
          },
          "grid"_a,
          "order"_a,
          "dimwise_global_mapping"_a = (r == 1) ? false : true,
          "dim_range"_a);

    return c;
  } // ... bind(...)
}; // class DiscontinuousLagrangeSpace


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_SPACES_L2_DISCONTINUOUS_LAGRANGE_HH
