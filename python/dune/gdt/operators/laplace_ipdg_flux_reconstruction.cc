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
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/operators/laplace-ipdg-flux-reconstruction.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/gdt/operators/interfaces.hh>


namespace Dune {
namespace GDT {
namespace bindings {


//template <class M, class MT, class GV>
//class LaplaceIpdgFluxReconstructionOperator
//{
//  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
//  static const size_t d = G::dimension;
//  using GP = XT::Grid::GridProvider<G>;

//public:
//  using type = GDT::LaplaceIpdgFluxReconstructionOperator<M, GV>;
//  using base_type = GDT::OperatorInterface<M, GV, 1, 1, d>;
//  using bound_type = pybind11::class_<type, base_type>;

//private:
//  using SS = typename type::SourceSpaceType;
//  using RS = typename type::RangeSpaceType;
//  using E = typename type::E;

//public:
//  static bound_type bind(pybind11::module& m,
//                         const std::string& matrix_id,
//                         const std::string& grid_id,
//                         const std::string& layer_id = "",
//                         const std::string& class_id = "laplace_ipdg_flux_reconstruction_operator")
//  {
//    namespace py = pybind11;
//    using namespace pybind11::literals;

//    const auto ClassName = XT::Common::to_camel_case(
//        bindings::OperatorInterface<M, GV, 1, d>::class_name(matrix_id, grid_id, layer_id, class_id));
//    bound_type c(m, ClassName.c_str(), ClassName.c_str());
//    c.def(py::init([](GP& grid,
//                      const SS& source_space,
//                      const RS& range_space,
//                      const double& symmetry_prefactor,
//                      const double& inner_penalty,
//                      const double& dirichlet_penalty,
//                      XT::Functions::GridFunction<E, d, d> diffusion,
//                      XT::Functions::GridFunction<E, d, d> weight_function) {
//            return new type(grid.leaf_view(),
//                            source_space,
//                            range_space,
//                            symmetry_prefactor,
//                            inner_penalty,
//                            dirichlet_penalty,
//                            diffusion,
//                            weight_function);
//          }),
//          "grid"_a,
//          "source_space"_a,
//          "range_space"_a,
//          "symmetry_prefactor"_a,
//          "inner_penalty"_a,
//          "dirichlet_penalty"_a,
//          "diffusion"_a,
//          "weight_function"_a,
//          py::keep_alive<1, 2>(),
//          py::keep_alive<1, 3>(),
//          py::keep_alive<1, 4>(),
//          py::keep_alive<1, 8>(),
//          py::keep_alive<1, 9>());

//    const auto FactoryName = XT::Common::to_camel_case(class_id);
//    if (std::is_same<MT, XT::LA::bindings::Istl>::value)
//      m.def(
//          FactoryName.c_str(),
//          [](GP& grid,
//             const SS& source_space,
//             const RS& range_space,
//             const double& symmetry_prefactor,
//             const double& inner_penalty,
//             const double& dirichlet_penalty,
//             XT::Functions::GridFunction<E, d, d> diffusion,
//             XT::Functions::GridFunction<E, d, d> weight_function,
//             const MT&) {
//            return new type(grid.leaf_view(),
//                            source_space,
//                            range_space,
//                            symmetry_prefactor,
//                            inner_penalty,
//                            dirichlet_penalty,
//                            diffusion,
//                            weight_function);
//          },
//          "grid"_a,
//          "source_space"_a,
//          "range_space"_a,
//          "symmetry_prefactor"_a,
//          "inner_penalty"_a,
//          "dirichlet_penalty"_a,
//          "diffusion"_a,
//          "weight_function"_a,
//          "la_backend"_a = MT(),
//          py::keep_alive<1, 2>(),
//          py::keep_alive<1, 3>(),
//          py::keep_alive<1, 4>(),
//          py::keep_alive<1, 8>(),
//          py::keep_alive<1, 9>());
//    else
//      m.def(
//          FactoryName.c_str(),
//          [](GP& grid,
//             const SS& source_space,
//             const RS& range_space,
//             const double& symmetry_prefactor,
//             const double& inner_penalty,
//             const double& dirichlet_penalty,
//             XT::Functions::GridFunction<E, d, d> diffusion,
//             XT::Functions::GridFunction<E, d, d> weight_function,
//             const MT&) {
//            return new type(grid.leaf_view(),
//                            source_space,
//                            range_space,
//                            symmetry_prefactor,
//                            inner_penalty,
//                            dirichlet_penalty,
//                            diffusion,
//                            weight_function);
//          },
//          "grid"_a,
//          "source_space"_a,
//          "range_space"_a,
//          "symmetry_prefactor"_a,
//          "inner_penalty"_a,
//          "dirichlet_penalty"_a,
//          "diffusion"_a,
//          "weight_function"_a,
//          "la_backend"_a,
//          py::keep_alive<1, 2>(),
//          py::keep_alive<1, 3>(),
//          py::keep_alive<1, 4>(),
//          py::keep_alive<1, 8>(),
//          py::keep_alive<1, 9>());

//    return c;
//  } // ... bind(...)
//}; // class LaplaceIpdgFluxReconstructionOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


//template <class M, class MT, class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
//struct LaplaceIpdgFluxReconstructionOperator_for_all_grids
//{
//  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
//  using GV = typename G::LeafGridView;

//  static void bind(pybind11::module& m, const std::string& matrix_id)
//  {
//    using Dune::GDT::bindings::LaplaceIpdgFluxReconstructionOperator;
//    using Dune::XT::Grid::bindings::grid_name;

//    LaplaceIpdgFluxReconstructionOperator<M, MT, GV>::bind(m, matrix_id, grid_name<G>::value());

//    LaplaceIpdgFluxReconstructionOperator_for_all_grids<M, MT, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(
//        m, matrix_id);
//  }
//};

//template <class M, class MT>
//struct LaplaceIpdgFluxReconstructionOperator_for_all_grids<M, MT, Dune::XT::Common::tuple_null_type>
//{
//  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
//};


PYBIND11_MODULE(_operators_laplace_ipdg_flux_reconstruction, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._operators_interfaces_common");
  py::module::import("dune.gdt._operators_interfaces_eigen");
  py::module::import("dune.gdt._operators_interfaces_istl_1d");
  py::module::import("dune.gdt._operators_interfaces_istl_2d");
  py::module::import("dune.gdt._operators_interfaces_istl_3d");
  py::module::import("dune.gdt._spaces_interface");

  /// \todo Add other linear algebra backends, if required!
//  LaplaceIpdgFluxReconstructionOperator_for_all_grids<LA::IstlRowMajorSparseMatrix<double>,
//                                                      LA::bindings::Istl,
//                                                      XT::Grid::bindings::AvailableGridTypes>::bind(m, "istl_sparse");
}
