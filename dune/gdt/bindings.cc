// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <config.h>

#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/string.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt.hh>

#include <dune/gdt/spaces.pbh>
#include <dune/gdt/assembler/system.pbh>
#include <dune/gdt/discretefunction.pbh>
#include <dune/gdt/functionals/l2.pbh>
#include <dune/gdt/operators/elliptic.pbh>
#include <dune/gdt/projections/dirichlet.pbh>

namespace py = pybind11;
using namespace pybind11::literals;


template <class SP>
void addbind_for_space(py::module& m,
                       const std::string& grid_id,
                       const std::string& layer_id,
                       const std::string& space_id,
                       const std::string& backend)
{
  using namespace Dune::XT;
  using Common::to_string;
  using Common::to_lower;

  typedef typename SP::type S;
  typedef typename Grid::extract_grid<typename S::GridViewType>::type G;
  typedef typename S::EntityType E;
  typedef typename S::DomainFieldType D;
  static const size_t d = S::dimDomain;
  typedef typename S::RangeFieldType R;
  typedef typename Dune::XT::LA::Container<R, Dune::XT::LA::Backends::istl_sparse>::MatrixType M;
  typedef typename Dune::XT::LA::Container<R, Dune::XT::LA::Backends::istl_sparse>::VectorType V;
  static const size_t r = S::dimRange;
  static const size_t rC = S::dimRangeCols;
  typedef Dune::XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1, 1> ScalarFunction;
  typedef Dune::XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunction;
  const std::string r_ = to_string(r);
  const std::string rC_ = to_string(rC);
  const std::string p_ = to_string(int(S::polOrder)); // without the int(...) we get linker errors on module import
  const std::string space_suffix = r_ + "x" + rC_ + "__p" + p_ + backend;
  // Space
  Dune::GDT::bind_space<S>(m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix);
  m.def(std::string("make_" + to_lower(space_id) + "_space__" + space_suffix + "__" + layer_id).c_str(),
        [](Grid::GridProvider<G>& grid_provider, const int level = 0) { return SP::create(grid_provider, level); },
        "grid_provider"_a,
        "level"_a = 0,
        py::keep_alive<0, 1>());
  // DiscreteFunction
  Dune::GDT::bind_DiscreteFunction<S, V>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
  // SystemAssembler
  Dune::GDT::bind_system_assembler<S>(m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix);
  // EllipticMatrixOperator
  Dune::GDT::bind_elliptic_matrix_operator<ScalarFunction, TensorFunction, S, M>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
  Dune::GDT::bind_elliptic_matrix_operator<ScalarFunction, void, S, M>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
  // L2VolumeVectorFunctional
  Dune::GDT::bind_l2_volume_vector_functional<ScalarFunction, S, V>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
  // L2FaceVectorFunctional
  Dune::GDT::bind_l2_face_vector_functional<ScalarFunction, S, V>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
} // ... addbind_for_space(...)


template <class SP>
void addbind_for_lagrange_space(py::module& m,
                                const std::string& grid_id,
                                const std::string& layer_id,
                                const std::string& space_id,
                                const std::string& backend)
{
  using namespace Dune::XT;
  using Common::to_string;

  typedef typename SP::type S;
  typedef typename S::EntityType E;
  typedef typename S::DomainFieldType D;
  static const size_t d = S::dimDomain;
  typedef typename S::RangeFieldType R;
  typedef typename Dune::XT::LA::Container<R, Dune::XT::LA::Backends::istl_sparse>::VectorType V;
  static const size_t r = S::dimRange;
  static const size_t rC = S::dimRangeCols;
  typedef Dune::XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1, 1> ScalarFunction;
  const std::string r_ = to_string(r);
  const std::string rC_ = to_string(rC);
  const std::string p_ = to_string(int(S::polOrder)); // without the int(...) we get linker errors on module import
  const std::string space_suffix = r_ + "x" + rC_ + "__p" + p_ + backend;

  // DirichletProjectionLocalizableOperator
  Dune::GDT::bind_DirichletProjectionLocalizableOperator<typename S::GridViewType,
                                                         ScalarFunction,
                                                         Dune::GDT::DiscreteFunction<S, V>>(
      m, space_id + "Space__" + grid_id + "_" + layer_id + "_to_" + space_suffix, "istl_sparse");
} // ... addbind_for_lagrange_space(...)


template <class G>
void addbind_for_grid(py::module& m, const std::string& grid_id)
{
  // FV
  addbind_for_space<Dune::GDT::FvSpaceProvider<G,
                                               Dune::XT::Grid::Layers::leaf,
                                               Dune::GDT::ChooseSpaceBackend::gdt,
                                               double,
                                               1,
                                               1>>(m, grid_id, "leaf", "Fv", "");
  addbind_for_space<Dune::GDT::FvSpaceProvider<G,
                                               Dune::XT::Grid::Layers::level,
                                               Dune::GDT::ChooseSpaceBackend::gdt,
                                               double,
                                               1,
                                               1>>(m, grid_id, "level", "Fv", "");
// CG
#if HAVE_DUNE_FEM
  addbind_for_space<Dune::GDT::CgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::leaf,
                                               Dune::GDT::ChooseSpaceBackend::fem,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "leaf", "Cg", "__fem");
  addbind_for_lagrange_space<Dune::GDT::CgSpaceProvider<G,
                                                        Dune::XT::Grid::Layers::leaf,
                                                        Dune::GDT::ChooseSpaceBackend::fem,
                                                        1,
                                                        double,
                                                        1,
                                                        1>>(m, grid_id, "leaf", "Cg", "__fem");
  addbind_for_space<Dune::GDT::CgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::level,
                                               Dune::GDT::ChooseSpaceBackend::fem,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "level", "Cg", "__fem");
  addbind_for_lagrange_space<Dune::GDT::CgSpaceProvider<G,
                                                        Dune::XT::Grid::Layers::level,
                                                        Dune::GDT::ChooseSpaceBackend::fem,
                                                        1,
                                                        double,
                                                        1,
                                                        1>>(m, grid_id, "level", "Cg", "__fem");
#endif
#if HAVE_DUNE_PDELAB
  addbind_for_space<Dune::GDT::CgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::leaf,
                                               Dune::GDT::ChooseSpaceBackend::pdelab,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "leaf", "Cg", "__pdelab");
  addbind_for_space<Dune::GDT::CgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::level,
                                               Dune::GDT::ChooseSpaceBackend::pdelab,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "level", "Cg", "__pdelab");
#endif
// DG
#if HAVE_DUNE_FEM
  addbind_for_space<Dune::GDT::DgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::leaf,
                                               Dune::GDT::ChooseSpaceBackend::fem,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "leaf", "Dg", "__fem");
  addbind_for_space<Dune::GDT::DgSpaceProvider<G,
                                               Dune::XT::Grid::Layers::level,
                                               Dune::GDT::ChooseSpaceBackend::fem,
                                               1,
                                               double,
                                               1,
                                               1>>(m, grid_id, "level", "Dg", "__fem");
#endif
  //#if HAVE_DUNE_PDELAB // these need to be guarded against simplicial grids
  //  addbind_for_space<Dune::GDT::DgSpaceProvider<G,
  //                                        Dune::XT::Grid::Layers::leaf,
  //                                        Dune::GDT::ChooseSpaceBackend::pdelab,
  //                                        1,
  //                                        double,
  //                                        1,
  //                                        1>>(m, grid_id, "leaf", "Dg", "__pdelab");
  //  addbind_for_space<Dune::GDT::DgSpaceProvider<G,
  //                                        Dune::XT::Grid::Layers::level,
  //                                        Dune::GDT::ChooseSpaceBackend::pdelab,
  //                                        1,
  //                                        double,
  //                                        1,
  //                                        1>>(m, grid_id, "level", "Dg", "__pdelab");
  //#endif
  // RT
  //#if HAVE_DUNE_PDELAB //  enabling those leads to problems in bindings of the functionals and vectors, those need to
  //  addbind_for_space<Dune::GDT::RtSpaceProvider<G, //                                 be disabled for bad dimensions
  //                                               Dune::XT::Grid::Layers::leaf,
  //                                               Dune::GDT::ChooseSpaceBackend::pdelab,
  //                                               0,
  //                                               double,
  //                                               G::dimension,
  //                                               1>>(m, grid_id, "leaf", "Rt", "__pdelab");
  //  addbind_for_space<Dune::GDT::RtSpaceProvider<G,
  //                                               Dune::XT::Grid::Layers::level,
  //                                               Dune::GDT::ChooseSpaceBackend::pdelab,
  //                                               0,
  //                                               double,
  //                                               G::dimension,
  //                                               1>>(m, grid_id, "level", "Rt", "__pdelab");
  //#endif
} // ... addbind_for_grid(...)


PYBIND11_PLUGIN(gdt)
{
  py::module m("gdt", "dune-gdt");

  py::module::import("common");
  py::module::import("grid");
  py::module::import("functions");
  py::module::import("la");

  m.def("init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = args.size();
          char** argv = Dune::XT::Common::vector_to_main_args(args);
#if HAVE_DUNE_FEM
          Dune::Fem::MPIManager::initialize(argc, argv);
#else
          Dune::MPIHelper::instance(argc, argv);
#endif
        },
        "args"_a = std::vector<std::string>());

  typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> G;
  const std::string grid_id = "2d_simplex_aluconform";

  addbind_for_grid<G>(m, grid_id);

  m.def("init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = -1,
        "max_debug_level"_a = -1,
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  return m.ptr();
}


#endif // HAVE_DUNE_PYBINDXI
