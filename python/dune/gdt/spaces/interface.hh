// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_SPACES_INTERFACE_HH
#define PYTHON_DUNE_GDT_SPACES_INTERFACE_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/gdt/spaces/interface.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, size_t r = 1, size_t rC = 1, class R = double>
class SpaceInterface
{
  using G = typename GV::Grid;
  static const size_t d = G::dimension;

public:
  using type = GDT::SpaceInterface<GV, r, rC, R>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "space_interface")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    if (r > 1)
      class_name += "_to_" + XT::Common::to_string(r) + "d";
    if (rC > 1)
      class_name += "x" + XT::Common::to_string(rC) + "d";
    if (!std::is_same<R, double>::value)
      class_name += "_" + XT::Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def_property_readonly("dimDomain", [](type&) { return d; });
    c.def_property_readonly("type", [](type& self) {
      std::stringstream ss;
      ss << self.type();
      return ss.str();
    });
    c.def_property_readonly("num_DoFs", [](type& self) { return self.mapper().size(); });
    c.def_property_readonly("min_polorder", [](type& self) { return self.min_polorder(); });
    c.def_property_readonly("max_polorder", [](type& self) { return self.max_polorder(); });
    c.def_property_readonly("is_lagrangian", [](type& self) { return self.is_lagrangian(); });
    c.def(
        "continuous", [](type& self, const int diff_order) { return self.continuous(diff_order); }, "diff_order"_a = 0);
    c.def("pre_adapt", [](type& self) { self.pre_adapt(); });
    c.def("adapt", [](type& self) { self.adapt(); });
    c.def("post_adapt", [](type& self) { return self.post_adapt(); });
    c.def("interpolation_points",
          [](type& self) {
            DUNE_THROW_IF(!self.is_lagrangian(), XT::Common::Exceptions::wrong_input_given, "");
            const auto& global_basis = self.basis();
            auto basis = global_basis.localize();
            DynamicVector<size_t> global_DoF_indices(self.mapper().max_local_size());
            auto points = std::make_unique<XT::LA::CommonDenseMatrix<double>>(self.mapper().size(), size_t(d), 0.);
            for (auto&& element : elements(self.grid_view())) {
              basis->bind(element);
              self.mapper().global_indices(element, global_DoF_indices);
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
    return c;
  } // ... bind(...)
}; // class SpaceInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_SPACES_INTERFACE_HH
