// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2023 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/xt/grid/dd/glued.hh>

using namespace Dune;

using G = YASP_2D_EQUIDISTANT_OFFSET;
using MG = YASP_2D_EQUIDISTANT_OFFSET;
using GridGlueType = XT::Grid::DD::Glued<MG, G, XT::Grid::Layers::leaf>;


PYBIND11_MODULE(_tools_interpolate_coarse_basis, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt");

  m.def("interpolate_coarse_basis",
        [](GridGlueType& dd_grid, const size_t I) {
          const auto& coarse_grid_view = dd_grid.macro_grid_view();
          // TODO: make coarse_space optional !
          auto coarse_space = GDT::make_discontinuous_lagrange_space(coarse_grid_view, /*order=*/1);
          auto coarse_basis = coarse_space.basis().localize();
          std::vector<XT::LA::IstlDenseVector<double>> interpolated_basis;
          for (auto&& macro_element : elements(coarse_grid_view)) {
            if (dd_grid.subdomain(macro_element) == I) {
              // this is the subdomain we are interested in, create space
               auto subdomain_grid_view = dd_grid.local_grid(macro_element).leaf_view();
               auto subdomain_space = GDT::make_continuous_lagrange_space(subdomain_grid_view, /*order=*/1);
              coarse_basis->bind(macro_element);
              for (size_t ii = 0; ii < coarse_basis->size(); ++ii)
                interpolated_basis.push_back(Dune::GDT::default_interpolation<XT::LA::IstlDenseVector<double>>(
                                                 coarse_basis->order(),
                                                 [&](const auto& point_in_physical_coordinates, const auto&) {
                                                   const auto point_macro_reference_element =
                                                       macro_element.geometry().local(point_in_physical_coordinates);
                                                   return coarse_basis->evaluate_set(point_macro_reference_element)[ii];
                                                 },
                                                 subdomain_space)
                                                 .dofs()
                                                 .vector());
              break;
            }
          }
          DUNE_THROW_IF(interpolated_basis.size() == 0, InvalidStateException, "This should not happen, I = " << I);
          return interpolated_basis;
        },
        py::call_guard<py::gil_scoped_release>(),
        "dd_grid"_a,
        "I"_a);
} // PYBIND11_MODULE(interpolate_coarse_basis, ...)
