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
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/operators/bilinear-form.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, size_t s_r = 1, size_t r_r = s_r>
class BilinearForm
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::BilinearForm<GV, s_r, 1, double, double, r_r>;
  using base_type = XT::Grid::ElementAndIntersectionFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  using E = typename type::E;
  using I = typename type::I;
  using F = typename type::ResultType;

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "bilinear_form")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r_r) + "d_range";
    class_name += "_" + XT::Common::to_string(s_r) + "d_source";
    const auto ClassName = XT::Common::to_camel_case(class_name);

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](GP& grid,
                      XT::Functions::GridFunction<E, s_r> source,
                      XT::Functions::GridFunction<E, r_r> range,
                      const std::string& logging_prefix) {
            return new type(grid.leaf_view(), source, range, logging_prefix);
          }),
          "grid"_a,
          "source"_a,
          "range"_a,
          "logging_prefix"_a = "",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>(),
          py::keep_alive<1, 4>());

    c.def_property_readonly("result", [](type& self) { return self.result(); });

    // methods
    c.def(
        "append",
        [](type& self,
           const LocalElementBilinearFormInterface<E, s_r, 1, F, F, r_r, 1, F>& local_bilinear_form,
           const XT::Common::Parameter& param,
           const XT::Grid::ElementFilter<GV>& filter) { self.append(local_bilinear_form, param, filter); },
        "local_element_bilinear_form"_a,
        "param"_a = XT::Common::Parameter(),
        "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type
           & (type::*)(const LocalElementBilinearFormInterface<E, s_r, 1, F, F, r_r, 1, F>&,
                       const XT::Common::Parameter&,
                       const XT::Grid::ElementFilter<GV>&))
              & type::append,
          "local_element_bilinear_form"_a,
          "param"_a = XT::Common::Parameter(),
          "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>(),
          py::is_operator());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type
           & (type::*)(const std::tuple<const LocalElementBilinearFormInterface<E, s_r, 1, F, F, r_r, 1, F>&,
                                        const XT::Common::Parameter&,
                                        const XT::Grid::ElementFilter<GV>&>&))
              & type::append,
          "tuple_of_localelementbilinearform_param_elementfilter"_a,
          py::is_operator());
    c.def(
        "apply2",
        [](type& self, const bool use_tbb) { return self.apply2(use_tbb); },
        "parallel"_a = false,
        py::call_guard<py::gil_scoped_release>());

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(
        FactoryName.c_str(),
        [](GP& grid,
           XT::Functions::GridFunction<E, s_r> source,
           XT::Functions::GridFunction<E, r_r> range,
           const std::string& logging_prefix) { return new type(grid.leaf_view(), source, range, logging_prefix); },
        "grid"_a,
        "source"_a,
        "range"_a,
        "logging_prefix"_a = "",
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class BilinearForm


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct BilinearForm_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::BilinearForm;
    using Dune::XT::Grid::bindings::grid_name;

    BilinearForm<GV>::bind(m, grid_name<G>::value());
    if (d > 1) {
      BilinearForm<GV, d, 1>::bind(m, grid_name<G>::value());
      BilinearForm<GV, 1, d>::bind(m, grid_name<G>::value());
      BilinearForm<GV, d, d>::bind(m, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...
    BilinearForm_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct BilinearForm_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_operators_bilinear_form, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_bilinear_forms_element_interface");

  BilinearForm_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
