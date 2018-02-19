// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.hh>

#include <python/dune/gdt/operators/l2.hh>
#include <dune/gdt/playground/spaces/restricted.hh>


#define DUNE_GDT_OPERATORS_L2_BIND(_m, _G, _gl, _s_backend, _s_type, _p, _r, _la)                                      \
  Dune::GDT::bindings::L2MatrixOperator<_G,                                                                            \
                                        Dune::XT::Grid::Layers::_gl,                                                   \
                                        Dune::GDT::SpaceType::_s_type,                                                 \
                                        Dune::GDT::Backends::_s_backend,                                               \
                                        _p,                                                                            \
                                        _r,                                                                            \
                                        Dune::XT::LA::Backends::_la>::bind(m)


using namespace Dune;
namespace py = pybind11;
using namespace pybind11::literals;
using Dune::XT::Grid::Layers;
using namespace Dune::XT;
using Dune::GDT::SpaceType;


template <class G,
          XT::Grid::Layers layer_type,
          XT::Grid::Backends layer_backend,
          size_t range_r = 1,
          size_t range_rC = 1,
          size_t source_r = range_r,
          size_t source_rC = range_rC>
void bind_l2_localizable_product(py::module& m)
{
  try {
    GDT::bindings::L2LocalizableProduct<G, layer_type, layer_backend, range_r, range_rC, source_r, source_rC>::bind(m);
  } catch (std::runtime_error&) {
  }
}


PYBIND11_MODULE(__operators_l2, m)
{

#if HAVE_DUNE_ALUGRID
  using G = ALU_2D_SIMPLEX_CONFORMING;

  bind_l2_localizable_product<G, Layers::dd_subdomain, XT::Grid::Backends::view>(m);
  bind_l2_localizable_product<G, Layers::leaf, XT::Grid::Backends::view>(m);
  DUNE_GDT_OPERATORS_L2_BIND(m, G, leaf, gdt, dg, 1, 1, istl_sparse);
  DUNE_GDT_OPERATORS_L2_BIND(m, G, leaf, gdt, dg, 2, 1, istl_sparse);
  DUNE_GDT_OPERATORS_L2_BIND(m, G, leaf, gdt, dg, 3, 1, istl_sparse);
  DUNE_GDT_OPERATORS_L2_BIND(m, G, dd_subdomain, gdt, block_dg, 1, 1, istl_sparse);
  DUNE_GDT_OPERATORS_L2_BIND(m, G, dd_subdomain, gdt, dg, 1, 1, istl_sparse);
  Dune::GDT::bindings::internal::
      L2MatrixOperator<GDT::RestrictedSpace<
                           typename GDT::
                               SpaceProvider<G, Layers::leaf, GDT::SpaceType::rt, GDT::Backends::gdt, 0, double, 2>::
                                   type,
                           typename XT::Grid::Layer<G,
                                                    Layers::dd_subdomain,
                                                    XT::Grid::Backends::view,
                                                    XT::Grid::DD::SubdomainGrid<G>>::type>,
                       XT::LA::IstlRowMajorSparseMatrix<double>>::bind(m,
                                                                       "RtAlu2dSimplexLeafRestrictedSubdomainPartSpace",
                                                                       "istl_row_major_sparse_matrix_double");

#endif // HAVE_DUNE_ALUGRID
  add_initialization(m, "dune.gdt.operators.elliptic");
}

#endif // HAVE_DUNE_PYBINDXI
