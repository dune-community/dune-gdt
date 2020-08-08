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

#include <dune/gdt/operators/localizable-operator.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>

#include "interfaces.hh"


namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class MT, class GV, size_t s_r = 1, size_t r_r = s_r>
class Operator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::LocalizableOperator<M, GV, s_r, 1, r_r>;
  using base_type = GDT::OperatorInterface<M, GV, s_r, 1, r_r>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using Eop = typename type::LocalElementOperatorType;

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        bindings::OperatorInterface<M, GV, s_r, r_r>::class_name(matrix_id, grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](const GP& grid, const SS& source_space, const RS& range_space, const bool parallel) {
            return new type(grid.leaf_view(), source_space, range_space, parallel);
          }),
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "parallel"_a = false,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>(),
          py::keep_alive<1, 4>());


    // methods
    c.def(
        "append",
        [](type& self, const Eop& local_op, const XT::Grid::ElementFilter<GV>& filter) {
          self.append(local_op, filter);
        },
        "local_element_operator"_a,
        "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type & (type::*)(const Eop&)) & type::operator+=,
          "local_element_operator"_a,
          py::is_operator());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type & (type::*)(const std::tuple<const Eop&, const XT::Grid::ElementFilter<GV>&>&)) & type::operator+=,
          "tuple_of_localelementop_elementfilter"_a,
          py::is_operator());
    /// \todo add intersection op

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    if (std::is_same<MT, XT::LA::bindings::Istl>::value)
      m.def(
          FactoryName.c_str(),
          [](const GP& grid, const SS& source_space, const RS& range_space, const bool parallel, const MT&) {
            return new type(grid.leaf_view(), source_space, range_space, parallel);
          },
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "parallel"_a = false,
          "la_backend"_a = MT(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());
    else
      m.def(
          FactoryName.c_str(),
          [](const GP& grid, const SS& source_space, const RS& range_space, const bool parallel, const MT&) {
            return new type(grid.leaf_view(), source_space, range_space, parallel);
          },
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "parallel"_a = false,
          "la_backend"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class Operator


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class M, class MT, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct Operator_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& matrix_id)
  {
    using Dune::GDT::bindings::Operator;
    using Dune::XT::Grid::bindings::grid_name;

    Operator<M, MT, GV>::bind(m, matrix_id, grid_name<G>::value());
    if (d > 1) {
      Operator<M, MT, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value());
      Operator<M, MT, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value());
      Operator<M, MT, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...
    Operator_for_all_grids<M, MT, typename GridTypes::tail_type>::bind(m, matrix_id);
  }
};

template <class M, class MT>
struct Operator_for_all_grids<M, MT, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
};


PYBIND11_MODULE(_operators_operator, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._local_operators_element_interface");
  py::module::import("dune.gdt._operators_interfaces_common");
  py::module::import("dune.gdt._operators_interfaces_eigen");
  py::module::import("dune.gdt._operators_interfaces_istl_1d");
  py::module::import("dune.gdt._operators_interfaces_istl_2d");
  py::module::import("dune.gdt._operators_interfaces_istl_3d");

  /// \todo Add other la backends if required
  Operator_for_all_grids<LA::IstlRowMajorSparseMatrix<double>, LA::bindings::Istl, XT::Grid::AvailableGridTypes>::bind(
      m, "istl_sparse");
}
