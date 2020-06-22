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
    Dune::GDT::bindings::OperatorInterface<M, G>::bind(m);
    Dune::GDT::bindings::ConstLincombOperator<M, G>::bind(m);
    Dune::GDT::bindings::LincombOperator<M, G>::bind(m);
    if (d > 1) {
      Dune::GDT::bindings::OperatorInterface<M, G, d, 1, 1, 1>::bind(m);
      Dune::GDT::bindings::ConstLincombOperator<M, G, d, 1, 1, 1>::bind(m);
      Dune::GDT::bindings::LincombOperator<M, G, d, 1, 1, 1>::bind(m);
      Dune::GDT::bindings::OperatorInterface<M, G, 1, 1, d, 1>::bind(m);
      Dune::GDT::bindings::ConstLincombOperator<M, G, 1, 1, d, 1>::bind(m);
      Dune::GDT::bindings::LincombOperator<M, G, 1, 1, d, 1>::bind(m);
      Dune::GDT::bindings::OperatorInterface<M, G, d, 1, d, 1>::bind(m);
      Dune::GDT::bindings::ConstLincombOperator<M, G, d, 1, d, 1>::bind(m);
      Dune::GDT::bindings::LincombOperator<M, G, d, 1, d, 1>::bind(m);
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
