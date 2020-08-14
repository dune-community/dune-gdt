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

#include <dune/gdt/local/integrands/laplace.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class G, class E, size_t r = 1, class F = double>
class LocalLaplaceIntegrand
{
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalLaplaceIntegrand<E, r, F>;
  using base_type = GDT::LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "local_laplace_integrand")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r) + "d_bases";
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<XT::Functions::GridFunction<E, d, d, F>, const std::string&>(),
          "diffusion"_a,
          "logging_prefix"_a = "",
          py::keep_alive<1, 2>());

    // factory
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    if (r == 1)
      m.def(
          FactoryName.c_str(),
          [](XT::Functions::GridFunction<E, d, d, F> diffusion,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(diffusion, logging_prefix); },
          "diffusion"_a,
          "dim_range_bases"_a = XT::Grid::bindings::Dimension<r>(),
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>());
    else
      m.def(
          FactoryName.c_str(),
          [](XT::Functions::GridFunction<E, d, d, F> diffusion,
             const XT::Grid::bindings::Dimension<r>&,
             const std::string& logging_prefix) { return new type(diffusion, logging_prefix); },
          "diffusion"_a,
          "dim_range_bases"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)
}; // class LocalLaplaceIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalLaplaceIntegrand_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::LocalLaplaceIntegrand;

    LocalLaplaceIntegrand<G, E>::bind(m);
    if (d > 1)
      LocalLaplaceIntegrand<G, E, d>::bind(m);
    // add your extra dimensions here
    // ...
    LocalLaplaceIntegrand_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalLaplaceIntegrand_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_laplace, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt._local_integrands_binary_element_interface");

  LocalLaplaceIntegrand_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
