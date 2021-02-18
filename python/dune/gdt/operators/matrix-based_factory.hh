// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_FACTORY_HH
#define DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_FACTORY_HH

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

#include "matrix-based.hh"

// The type itself is bound by the OperatorInterface, required as return type for jacobian
template <class M, class MT, class ST, class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct MatrixOperatorFactory_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& matrix_id)
  {
    using Dune::GDT::bindings::MatrixOperator;
    using Dune::XT::Grid::bindings::grid_name;

    MatrixOperator<M, MT, ST, LGV>::bind_leaf_factory(m, matrix_id);
    if (d > 1) {
      MatrixOperator<M, MT, ST, LGV, d, 1>::bind_leaf_factory(m, matrix_id);
      MatrixOperator<M, MT, ST, LGV, 1, d>::bind_leaf_factory(m, matrix_id);
      MatrixOperator<M, MT, ST, LGV, d, d>::bind_leaf_factory(m, matrix_id);
    }
    // add your extra dimensions here
    // ...
#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      MatrixOperator<M, MT, ST, CGV, 1, 1, LGV, LGV>::bind_coupling_factory(m, matrix_id);
      if (d > 1) {
        MatrixOperator<M, MT, ST, CGV, d, 1, LGV, LGV>::bind_coupling_factory(m, matrix_id);
        MatrixOperator<M, MT, ST, CGV, 1, d, LGV, LGV>::bind_coupling_factory(m, matrix_id);
        MatrixOperator<M, MT, ST, CGV, d, d, LGV, LGV>::bind_coupling_factory(m, matrix_id);
      }
    }
#endif
    MatrixOperatorFactory_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m, matrix_id);
  }
};

template <class M, class MT, class ST>
struct MatrixOperatorFactory_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
};

#endif // DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_FACTORY_HH
