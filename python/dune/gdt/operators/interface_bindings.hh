// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_OPERATORS_INTERFACE_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_INTERFACE_BINDINGS_HH

#include "interfaces.hh"
#include "lincomb.hh"


template <class M, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct OperatorInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::ConstLincombOperator;
    using Dune::GDT::bindings::LincombOperator;
    using Dune::GDT::bindings::OperatorInterface;
    using Dune::XT::Grid::bindings::grid_name;

    OperatorInterface<M, G>::bind(m, grid_name<G>::value());
    ConstLincombOperator<M, G>::bind(m, grid_name<G>::value());
    LincombOperator<M, G>::bind(m, grid_name<G>::value());
    if (d > 1) {
      OperatorInterface<M, G, d, 1, 1, 1>::bind(m, grid_name<G>::value());
      ConstLincombOperator<M, G, d, 1, 1, 1>::bind(m, grid_name<G>::value());
      LincombOperator<M, G, d, 1, 1, 1>::bind(m, grid_name<G>::value());
      OperatorInterface<M, G, 1, 1, d, 1>::bind(m, grid_name<G>::value());
      ConstLincombOperator<M, G, 1, 1, d, 1>::bind(m, grid_name<G>::value());
      LincombOperator<M, G, 1, 1, d, 1>::bind(m, grid_name<G>::value());
      OperatorInterface<M, G, d, 1, d, 1>::bind(m, grid_name<G>::value());
      ConstLincombOperator<M, G, d, 1, d, 1>::bind(m, grid_name<G>::value());
      LincombOperator<M, G, d, 1, d, 1>::bind(m, grid_name<G>::value());
    }
    // add your extra dimensions here
    // ...

    OperatorInterface_for_all_grids<M, typename GridTypes::tail_type>::bind(m);
  }
};

template <class M>
struct OperatorInterface_for_all_grids<M, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


#endif // PYTHON_DUNE_GDT_OPERATORS_INTERFACE_BINDINGS_HH
