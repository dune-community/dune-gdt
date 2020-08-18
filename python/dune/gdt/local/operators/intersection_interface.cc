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
#include <dune/gdt/local/operators/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/la/container.bindings.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class GV, size_t s_r, size_t r_r>
class LocalIntersectionOperatorInterface
{
  static_assert(XT::LA::is_vector<V>::value, "");
  static_assert(XT::Grid::is_view<GV>::value, "");

  using I = XT::Grid::extract_intersection_t<GV>;

public:
  using type = GDT::LocalIntersectionOperatorInterface<I, V, GV, s_r, 1, double, r_r, 1, double>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& vector_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "local_intersection_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += vector_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r_r) + "d_range";
    class_name += "_" + XT::Common::to_string(s_r) + "d_source";
    class_name += "_interface";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    // static information about dims
    // ...
    c.def("copy", [](const type& self) { return self.copy(); });
    c.def_property_readonly("linear", [](const type& self) { return self.linear(); });
    return c;
  } // ... bind(...)
}; // class LocalIntersectionOperatorInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalIntersectionOperatorInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& vector_id)
  {
    using Dune::GDT::bindings::LocalIntersectionOperatorInterface;
    using Dune::XT::Grid::bindings::grid_name;

    LocalIntersectionOperatorInterface<V, GV, 1, 1>::bind(m, vector_id, grid_name<G>::value());
    if (d > 1) {
      LocalIntersectionOperatorInterface<V, GV, 1, d>::bind(m, vector_id, grid_name<G>::value());
      LocalIntersectionOperatorInterface<V, GV, d, 1>::bind(m, vector_id, grid_name<G>::value());
      LocalIntersectionOperatorInterface<V, GV, d, d>::bind(m, vector_id, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...
    LocalIntersectionOperatorInterface_for_all_grids<V, typename GridTypes::tail_type>::bind(m, vector_id);
  }
};

template <class V>
struct LocalIntersectionOperatorInterface_for_all_grids<V, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*vector_id*/) {}
};


PYBIND11_MODULE(_local_operators_intersection_interface, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  LocalIntersectionOperatorInterface_for_all_grids<XT::LA::IstlDenseVector<double>, XT::Grid::AvailableGridTypes>::bind(
      m, "istl_dense");
}
