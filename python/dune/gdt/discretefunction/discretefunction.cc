// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include "config.h"

#include <dune/xt/grid/grids.hh>
#include <python/dune/xt/la/traits.hh>

#include "discretefunction.hh"


template <class V, class VT, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct DiscreteFunction_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::DiscreteFunction<V, VT, GV>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::DiscreteFunction<V, VT, GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    DiscreteFunction_for_all_grids<V, VT, typename GridTypes::tail_type>::bind(m);
  }
};

template <class V, class VT>
struct DiscreteFunction_for_all_grids<V, VT, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_discretefunction_discretefunction, m)
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
  py::module::import("dune.gdt._discretefunction_dof_vector");

  DiscreteFunction_for_all_grids<LA::CommonDenseVector<double>, LA::bindings::Common, XT::Grid::AvailableGridTypes>::
      bind(m);
#if HAVE_EIGEN
  DiscreteFunction_for_all_grids<LA::EigenDenseVector<double>, LA::bindings::Eigen, XT::Grid::AvailableGridTypes>::bind(
      m);
#endif
  DiscreteFunction_for_all_grids<LA::IstlDenseVector<double>, LA::bindings::Istl, XT::Grid::AvailableGridTypes>::bind(
      m);
}