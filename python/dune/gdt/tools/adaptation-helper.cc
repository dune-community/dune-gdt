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
#include <dune/xt/la/container/common/vector.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/adaptation-helper.hh>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>
#include <python/dune/xt/la/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


/**
 * \note In addition to providing the bindings, this is also a trampoline class to store the additional vector of
 *       markings
 */
template <class V, class VT, class GV, size_t r = 1, size_t rC = 1, class RF = double>
class AdaptationHelper : public GDT::AdaptationHelper<V, GV, r, rC, RF>
{
  using G = XT::Grid::extract_grid_t<GV>;
  static_assert(rC == 1, "");

  using ThisType = GDT::bindings::AdaptationHelper<V, VT, GV, r, rC, RF>;
  using BaseType = GDT::AdaptationHelper<V, GV, r, rC, RF>;

  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::SpaceType;

public:
  using type = ThisType;
  using bound_type = pybind11::class_<type>;

  AdaptationHelper(G& grd, const std::string& logging_prefix = "")
    : BaseType(grd, logging_prefix)
    , marker_indices(grd.leafGridView()) // BAD, bad, bad, has to coincide with GV!
    , markers(marker_indices.mapper().size())
  {}

  GDT::FiniteVolumeSpace<GV, r, rC> marker_indices;
  XT::LA::CommonDenseVector<size_t> markers;

  G& grid() // expose access for bindings below
  {
    return this->grid_;
  }

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "adaptation_helper")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    if (r > 1)
      class_name += "_to_" + XT::Common::to_string(r) + "d";
    if (!std::is_same<RF, double>::value)
      class_name += "_" + XT::Common::Typename<RF>::value(/*fail_wo_typeid=*/true);
    class_name += "_" + XT::LA::bindings::container_name<V>::value();
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](XT::Grid::GridProvider<G>& grid, const std::string& logging_prefix) {
            return new type(grid.grid(), logging_prefix);
          }),
          "grid"_a,
          "logging_prefix"_a = "");

    c.def_readwrite("markers", &type::markers);

    c.def("append", [](type& self, SpaceType& space, DiscreteFunctionType& discrete_function) {
      self.append(space, discrete_function);
    });
    c.def(
        "mark",
        [](type& self) {
          auto& grid = self.grid();
          auto grid_view = grid.leafGridView();
          for (auto&& element : elements(grid_view))
            grid.mark(XT::Common::numeric_cast<int>(
                          self.markers[self.marker_indices.mapper().global_index(element, size_t(0))]),
                      element);
        },
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "pre_adapt",
        [](type& self, const bool pre_adapt_grid) { self.pre_adapt(pre_adapt_grid); },
        "pre_adapt_grid"_a = true,
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "adapt",
        [](type& self, const bool adapt_grid) { self.adapt(adapt_grid); },
        "adapt_grid"_a = true,
        py::call_guard<py::gil_scoped_release>());
    c.def(
        "post_adapt",
        [](type& self, const bool post_adapt_grid, const bool clear, const bool indicate_new_elements) {
          // we need to keep track of these elements before ...
          self.marker_indices.update_after_adapt();
          self.markers.resize(self.marker_indices.mapper().size());
          self.markers *= 0.;
          if (indicate_new_elements) {
            auto grid_view = self.grid().leafGridView();
            for (auto&& element : elements(grid_view))
              if (element.isNew())
                self.markers[self.marker_indices.mapper().global_index(element, size_t(0))] = size_t(1);
          }
          // ... possibly cleaning up the grid
          self.post_adapt(post_adapt_grid, clear);
        },
        "post_adapt_grid"_a = true,
        "clear"_a = false,
        "indicate_new_elements"_a = false,
        py::call_guard<py::gil_scoped_release>());

    const auto FactoryName = XT::Common::to_camel_case(class_id);
    if (std::is_same<VT, XT::LA::bindings::Istl>::value && r == 1)
      m.def(
          FactoryName.c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const VT&,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.grid(), logging_prefix); },
          "grid"_a,
          "la_backend"_a = XT::LA::bindings::Istl(),
          "dim_range"_a = XT::Grid::bindings::Dimension<r>(),
          "logging_prefix"_a = "");
    else if (std::is_same<VT, XT::LA::bindings::Istl>::value)
      m.def(
          FactoryName.c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const VT&,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.grid(), logging_prefix); },
          "grid"_a,
          "dim_range"_a,
          "la_backend"_a = XT::LA::bindings::Istl(),
          "logging_prefix"_a = "");
    else if (r == 1)
      m.def(
          FactoryName.c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const VT&,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.grid(), logging_prefix); },
          "grid"_a,
          "la_backend"_a,
          "dim_range"_a = XT::Grid::bindings::Dimension<r>(),
          "logging_prefix"_a = "");
    else
      m.def(
          FactoryName.c_str(),
          [](XT::Grid::GridProvider<G>& grid,
             const VT&,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(grid.grid(), logging_prefix); },
          "grid"_a,
          "la_backend"_a,
          "dim_range"_a,
          "logging_prefix"_a = "");

    return c;
  } // ... bind(...)
}; // class AdaptationHelper


} // namespace bindings
} // namespace GDT
} // namespace Dune


using AvailableAdaptiveGridTypes = std::tuple<ONED_1D,
#if HAVE_DUNE_ALUGRID
                                              ALU_2D_SIMPLEX_CONFORMING,
                                              ALU_3D_SIMPLEX_CONFORMING
#endif
                                              >;


template <class V, class VT, class GridTypes = AvailableAdaptiveGridTypes>
struct AdaptationHelper_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::AdaptationHelper;

    AdaptationHelper<V, VT, GV>::bind(m);
    if (d > 1)
      AdaptationHelper<V, VT, GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    AdaptationHelper_for_all_grids<V, VT, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <class V, class VT>
struct AdaptationHelper_for_all_grids<V, VT, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_tools_adaptation_helper, m)
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
  py::module::import("dune.gdt._discretefunction_discretefunction");

  //  AdaptationHelper_for_all_grids<LA::CommonDenseVector<double>, LA::bindings::Common>::bind(m);
  //#if HAVE_EIGEN
  //  AdaptationHelper_for_all_grids<LA::EigenDenseVector<double>, LA::bindings::Eigen>::bind(m);
  //#endif
  AdaptationHelper_for_all_grids<LA::IstlDenseVector<double>, LA::bindings::Istl>::bind(m);
}
