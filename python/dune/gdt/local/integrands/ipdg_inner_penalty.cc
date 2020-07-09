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

#include <dune/gdt/local/integrands/ipdg.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class I>
class LocalIPDGInnerPenaltyIntegrand
{
  using G = XT::Grid::extract_grid_t<I>;
  using E = XT::Grid::extract_entity_t<G>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;
  static const size_t r = 1;

public:
  using type = GDT::LocalIPDGIntegrands::InnerPenalty<I>;
  using base_type = GDT::LocalQuaternaryIntersectionIntegrandInterface<I>;
  using bound_type = pybind11::class_<type, base_type>;
  using F = typename type::F;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_IPDG_inner_penalty_integrand",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
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
    c.def(py::init<const double&, XT::Functions::GridFunction<E, d, d, F>>(),
          "penalty"_a,
          "weight"_a,
          py::keep_alive<1, 3>());

    // factory
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(FactoryName.c_str(),
          [](const double& penalty, XT::Functions::GridFunction<E, d, d, F> weight) { return type(penalty, weight); },
          "penalty"_a,
          "weight"_a,
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class LocalIPDGInnerPenaltyIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalIPDGInnerPenaltyIntegrand_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalIPDGInnerPenaltyIntegrand<I>::bind(m);

    LocalIPDGInnerPenaltyIntegrand_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct LocalIPDGInnerPenaltyIntegrand_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_ipdg_inner_penalty, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt._local_integrands_quaternary_intersection_interface");

  LocalIPDGInnerPenaltyIntegrand_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}