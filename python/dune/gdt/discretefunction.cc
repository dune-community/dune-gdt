// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#  include <dune/pybindxi/pybind11.h>

#  include <python/dune/xt/common/exceptions.bindings.hh>
#  include <dune/xt/la/container.hh>
#  include <dune/xt/grid/grids.hh>

#  include "discretefunction.hh"


template <class V, class Tuple = Dune::XT::Grid::AvailableGridTypes>
struct binder_for_all_grids
{
  using G = typename Tuple::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::DiscreteFunction<V, GV>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::DiscreteFunction<V, GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    binder_for_all_grids<V, typename Tuple::tail_type>::bind(m);
  }
};

template <class V>
struct binder_for_all_grids<V, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class Tuple = Dune::XT::LA::AvailableVectorTypes<double>>
void bind_for_all_vectors_and_grids(pybind11::module& m)
{
  using V = typename Tuple::head_type;
  binder_for_all_grids<V>::bind(m);
  bind_for_all_vectors_and_grids<typename Tuple::tail_type>(m);
}

template <>
void bind_for_all_vectors_and_grids<boost::tuples::null_type>(pybind11::module& /*m*/)
{}


PYBIND11_MODULE(discretefunction, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  Dune::XT::Common::bindings::addbind_exceptions(m);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt");

  bind_for_all_vectors_and_grids<>(m);
} // PYBIND11_MODULE(discretefunction, ...)


#endif // HAVE_DUNE_PYBINDXI
