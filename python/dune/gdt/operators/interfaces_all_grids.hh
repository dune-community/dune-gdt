// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH

#include "interfaces.hh"
#include "lincomb.hh"
#include "matrix-based.hh"


/**
 * We require ConstLincombOperator and LincombOperator for the numeric operators and MatrixOperator for the jacobian.
 */
template <class M, class MT, class ST, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct OperatorInterface_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& matrix_id)
  {
    using Dune::GDT::bindings::ConstLincombOperator;
    using Dune::GDT::bindings::LincombOperator;
    using Dune::GDT::bindings::MatrixOperator;
    using Dune::GDT::bindings::OperatorInterface;
    using Dune::XT::Grid::bindings::grid_name;

    OperatorInterface<M, GV>::bind(m, matrix_id, grid_name<G>::value());
    ConstLincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value());
    LincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value());
    MatrixOperator<M, MT, ST, GV>::bind_type(m, matrix_id, grid_name<G>::value());
    if (d > 1) {
      OperatorInterface<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value());
      ConstLincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value());
      LincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value());
      MatrixOperator<M, MT, ST, GV, d, 1>::bind_type(m, matrix_id, grid_name<G>::value());

      OperatorInterface<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value());
      ConstLincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value());
      LincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value());
      MatrixOperator<M, MT, ST, GV, 1, d>::bind_type(m, matrix_id, grid_name<G>::value());

      OperatorInterface<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
      ConstLincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
      LincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
      MatrixOperator<M, MT, ST, GV, d, d>::bind_type(m, matrix_id, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...

    OperatorInterface_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m, matrix_id);
  }
};

template <class M, class MT, class ST>
struct OperatorInterface_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
};


#endif // PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH
