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
#include <dune/gdt/local/operators/indicator.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/la/container.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class GV, size_t s_r = 1>
class LocalIntersectionBilinearFormIndicatorOperator
{
  static_assert(XT::LA::is_vector<V>::value, "");
  static_assert(XT::Grid::is_view<GV>::value, "");
  using I = XT::Grid::extract_intersection_t<GV>;

public:
  using type = GDT::LocalIntersectionBilinearFormIndicatorOperator<I, V, GV, s_r>;
  using base_type = typename type::BaseType;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& vector_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "local_intersection_bilinear_form_indicator_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += vector_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(s_r) + "d_source";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const typename type::LocalBilinearFormType&>(), "local_intersection_bilinear_form"_a);
    c.def(py::init<const typename type::LocalBilinearFormType&, const typename type::SourceType&>(),
          "local_intersection_bilinear_form"_a,
          "source"_a);

    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const typename type::LocalBilinearFormType& local_intersection_bilinear_form,
           const typename type::SourceType& source) { return new type(local_intersection_bilinear_form, source); },
        "local_intersection_bilinear_form"_a,
        "source"_a);
    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const typename type::LocalBilinearFormType& local_intersection_bilinear_form) {
          return new type(local_intersection_bilinear_form);
        },
        "local_intersection_bilinear_form"_a);

    return c;
  } // ... bind(...)
}; // class LocalIntersectionBilinearFormIndicatorOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalIntersectionBilinearFormIndicatorOperator_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& vector_id)
  {
    using Dune::GDT::bindings::LocalIntersectionBilinearFormIndicatorOperator;
    using Dune::XT::Grid::bindings::grid_name;

    LocalIntersectionBilinearFormIndicatorOperator<V, GV>::bind(m, vector_id, grid_name<G>::value());
    if (d > 1) {
      LocalIntersectionBilinearFormIndicatorOperator<V, GV, d>::bind(m, vector_id, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...
    LocalIntersectionBilinearFormIndicatorOperator_for_all_grids<V, typename GridTypes::tail_type>::bind(m, vector_id);
  }
};

template <class V>
struct LocalIntersectionBilinearFormIndicatorOperator_for_all_grids<V, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*vector_id*/) {}
};


PYBIND11_MODULE(_local_operators_intersection_indicator, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_bilinear_forms_intersection_interface");
  py::module::import("dune.gdt._local_operators_intersection_interface");

  LocalIntersectionBilinearFormIndicatorOperator_for_all_grids<XT::LA::IstlDenseVector<double>,
                                                               XT::Grid::AvailableGridTypes>::bind(m, "istl_dense");
}
