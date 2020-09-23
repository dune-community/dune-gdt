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

#include <dune/gdt/local/integrands/laplace-ipdg.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class G, class I>
class LocalLaplaceIPDGDirichletCouplingIntegrand
{
  using E = XT::Grid::extract_entity_t<G>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalLaplaceIPDGIntegrands::DirichletCoupling<I>;
  using base_unary_type = GDT::LocalUnaryIntersectionIntegrandInterface<I>;
  using base_binary_type = GDT::LocalBinaryIntersectionIntegrandInterface<I>;
  using bound_type = pybind11::class_<type, base_unary_type, base_binary_type>;
  using F = typename type::F;

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "local_laplace_IPDG_dirichlet_coupling_integrand")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_1d_bases";
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](const double& symmetry_prefactor,
                      const XT::Functions::GridFunctionInterface<E, d, d, F>& diffusion,
                      const std::string& logging_prefix) {
            return new type(symmetry_prefactor, diffusion, /*dirichlet_data=*/0., logging_prefix);
          }),
          "symmetry_prefactor"_a,
          "diffusion"_a,
          "logging_prefix"_a = "");
    c.def(py::init([](const double& symmetry_prefactor,
                      const XT::Functions::GridFunctionInterface<E, d, d, F>& diffusion,
                      const XT::Functions::GridFunctionInterface<E, 1, 1, F>& dirichlet_data,
                      const std::string& logging_prefix) {
            return new type(symmetry_prefactor, diffusion, dirichlet_data, logging_prefix);
          }),
          "symmetry_prefactor"_a,
          "diffusion"_a,
          "dirichlet_data"_a,
          "logging_prefix"_a = "");

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(
        FactoryName.c_str(),
        [](const double& symmetry_prefactor,
           const XT::Functions::GridFunctionInterface<E, d, d, F>& diffusion,
           const std::string& logging_prefix) {
          return new type(symmetry_prefactor, diffusion, /*dirichlet_data=*/0., logging_prefix);
        },
        "symmetry_prefactor"_a,
        "diffusion"_a,
        "logging_prefix"_a = "");
    m.def(
        FactoryName.c_str(),
        [](const double& symmetry_prefactor,
           const XT::Functions::GridFunctionInterface<E, d, d, F>& diffusion,
           const XT::Functions::GridFunctionInterface<E, 1, 1, F>& dirichlet_data,
           const std::string& logging_prefix) {
          return new type(symmetry_prefactor, diffusion, dirichlet_data, logging_prefix);
        },
        "symmetry_prefactor"_a,
        "diffusion"_a,
        "dirichlet_data"_a,
        "logging_prefix"_a = "");

    return c;
  } // ... bind(...)
}; // class LocalLaplaceIPDGDirichletCouplingIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct LocalLaplaceIPDGDirichletCouplingIntegrand_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::LocalLaplaceIPDGDirichletCouplingIntegrand<G, I>::bind(m);

    LocalLaplaceIPDGDirichletCouplingIntegrand_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct LocalLaplaceIPDGDirichletCouplingIntegrand_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_laplace_ipdg_dirichlet_coupling, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;
  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt._local_integrands_unary_intersection_interface");
  py::module::import("dune.gdt._local_integrands_binary_intersection_interface");

  LocalLaplaceIPDGDirichletCouplingIntegrand_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
