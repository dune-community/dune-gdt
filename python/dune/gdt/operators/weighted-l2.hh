// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_OPERATORS_WEIGHTED_L2_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_WEIGHTED_L2_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <python/dune/xt/grid/layers.bindings.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/operators/weighted-l2.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class G, XT::Grid::Layers layer_type, XT::Grid::Backends layer_backend>
class WeightedL2LocalizableProduct
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename XT::Grid::Layer<G, layer_type, layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;

  typedef XT::Grid::extract_entity_t<GL> E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, double, 1, 1> F;

  template <bool is_dd_subdomain = (layer_type == XT::Grid::Layers::dd_subdomain)
                                   || (layer_type == XT::Grid::Layers::dd_subdomain_boundary)
                                   || (layer_type == XT::Grid::Layers::dd_subdomain_coupling),
            bool anything = true>
  struct helper
  {
    static void bind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(std::string("apply_weighted_l2_product_" + XT::Grid::layer_names[layer_type] + "_"
                        + XT::Grid::bindings::backend_name<layer_backend>::value())
                .c_str(),
            [](const F& weight,
               const F& range,
               const F& source,
               const XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid,
               const int level_or_subdomain,
               const size_t over_integrate) {
              return GDT::WeightedL2LocalizableProduct<F, GL, F, F>(
                         over_integrate,
                         weight,
                         grid.template layer<layer_type, layer_backend>(level_or_subdomain),
                         range,
                         source)
                  .apply2();
            },
            "weight"_a,
            "range"_a,
            "source"_a,
            "grid"_a,
            "level_or_subdomain"_a = -1,
            "over_integrate"_a = 0);
    } // ... bind(...)
  }; // struct helper<true, ...>

  template <bool anything>
  struct helper<false, anything>
  {
    static void bind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(std::string("apply_weighted_l2_product_" + XT::Grid::layer_names[layer_type] + "_"
                        + XT::Grid::bindings::backend_name<layer_backend>::value())
                .c_str(),
            [](const F& weight,
               const F& range,
               const F& source,
               const XT::Grid::GridProvider<G, Dune::XT::Grid::none_t>& grid,
               const int level,
               const size_t over_integrate) {
              return GDT::WeightedL2LocalizableProduct<F, GL, F, F>(
                         over_integrate, weight, grid.template layer<layer_type, layer_backend>(level), range, source)
                  .apply2();
            },
            "weight"_a,
            "range"_a,
            "source"_a,
            "grid"_a,
            "level"_a = -1,
            "over_integrate"_a = 0);
      m.def(std::string("apply_weighted_l2_product_" + XT::Grid::layer_names[layer_type] + "_"
                        + XT::Grid::bindings::backend_name<layer_backend>::value())
                .c_str(),
            [](const F& weight,
               const F& range,
               const F& source,
               const XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid,
               const int layer,
               const size_t over_integrate) {
              return GDT::WeightedL2LocalizableProduct<F, GL, F, F>(
                         over_integrate, weight, grid.template layer<layer_type, layer_backend>(layer), range, source)
                  .apply2();
            },
            "weight"_a,
            "range"_a,
            "source"_a,
            "grid"_a,
            "layer"_a = -1,
            "over_integrate"_a = 0);
    } // ... bind(...)
  }; // struct helper<false, ...>

public:
  static void bind(pybind11::module& m)
  {
    helper<>::bind(m);
  }
}; // class WeightedL2LocalizableProduct


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_OPERATORS_WEIGHTED_L2_BINDINGS_HH
