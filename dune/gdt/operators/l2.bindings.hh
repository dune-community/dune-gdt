// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_L2_BINDINGS_HH
#define DUNE_GDT_OPERATORS_L2_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.hh>
#include <dune/gdt/spaces.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "base.bindings.hh"
#include "l2.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class G,
          XT::Grid::Layers layer_type,
          GDT::SpaceType space_type,
          GDT::Backends space_backend,
          int p,
          size_t r,
          XT::LA::Backends la_backend>
class L2MatrixOperator
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef GDT::SpaceProvider<G, layer_type, space_type, space_backend, p, double, r, 1> RP;
  typedef typename RP::type R;
  typedef typename XT::LA::Container<double, la_backend>::MatrixType M;

public:
  typedef GDT::L2MatrixOperator<R, M /*, GL, S, F*/> type;
  typedef pybind11::class_<type, GDT::SystemAssembler<R>> bound_type;

public:
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("l2_matrix_operator_" + space_name<RP>::value() + "_"
                                                     + XT::LA::bindings::container_name<M>::value());

    auto c = MatrixOperatorBase<type>::bind(m, ClassName);

    m.def(std::string("make_l2_matrix_operator_" + XT::LA::bindings::container_name<M>::value()).c_str(),
          [](const R& space, const size_t over_integrate) {
            return make_l2_matrix_operator<M>(space, over_integrate).release();
          },
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>());

    m.def("make_l2_matrix_operator",
          [](M& matrix, const R& space, const size_t over_integrate) {
            return make_l2_matrix_operator(matrix, space, over_integrate).release();
          },
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class L2MatrixOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_L2_BINDINGS_HH
