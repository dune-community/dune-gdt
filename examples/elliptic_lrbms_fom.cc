// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2022)

#include "config.h"

#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

//#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/grid-function.hh>

//#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
//#include <dune/gdt/local/functionals/integrals.hh>
//#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/laplace.hh>
//#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/bilinear-form.hh>
//#include <dune/gdt/operators/matrix.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>
#include <dune/gdt/test/stationary-heat-equation/ESV2007.hh>


using namespace Dune;
using namespace Dune::GDT;


// some global defines
// - macro grid
using G_M = YASP_2D_EQUIDISTANT_OFFSET;
// static const constexpr size_t d = G_M::dimension;
using GV_M = typename G_M::LeafGridView;
// - local subdomain grids (we use different grids on purpose, if available)
#if HAVE_DUNE_ALUGRID
using G_L = ALU_2D_SIMPLEX_CONFORMING;
#else
using G_L = YASP_2D_EQUIDISTANT_OFFSET;
#endif
using GV_L = typename G_L::LeafGridView;
using E_L = XT::Grid::extract_entity_t<GV_L>;
using I_L = XT::Grid::extract_intersection_t<GV_L>;
// - linear algebra
using M = XT::LA::IstlRowMajorSparseMatrix<double>;
using V = XT::LA::IstlDenseVector<double>;


int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    logger.info() << "Creating problem, DD grid, and local spaces ... " << std::flush;

    // analytical problem
    const double diffusion = 1;
    const double force = 1;

    // grid
    auto macro_grid = XT::Grid::make_cube_grid<G_M>(0., 1., 8);
    XT::Grid::DD::Glued<G_M, G_L> dd_grid(
        macro_grid,
        /*num_local_refinements=*/2,
        /*prepare_glues=*/false,
        /*allow_for_broken_orientation_of_coupling_intersections=*/true); // <- get fixed in coupling view anyway
    auto macro_grid_view = dd_grid.macro_grid_view();

    // local spaces
    std::vector<std::unique_ptr<SpaceInterface<GV_L>>> local_spaces;
    local_spaces.reserve(dd_grid.num_subdomains());
    for (size_t ss = 0; ss < dd_grid.num_subdomains(); ++ss)
      if (DXTC_CONFIG_GET("locally_conforming", false))
        local_spaces.emplace_back(
            new ContinuousLagrangeSpace<GV_L>(dd_grid.local_grid(ss).leaf_view(), DXTC_CONFIG_GET("local_order", 1)));
      else
        local_spaces.emplace_back(new DiscontinuousLagrangeSpace<GV_L>(dd_grid.local_grid(ss).leaf_view(),
                                                                       DXTC_CONFIG_GET("local_order", 1)));

    logger.info() << "done" << std::endl;
    logger.info() << "Creating global DG DoF mapping and global containers ... " << std::flush;

    // macro DG discretization
    const double symmetry_prefactor = 1; // SIPDG
    const double penalty_parameter = 16; // non-degenerate simplicial grids in 2d
    const auto& weight = diffusion; // SWIPDG, not SIPDG
    // For SWIPDG, the weight should be
    // - the diffusion; for non-parametric problems
    // - the diffusion for the parameter inducing the energy norm; for parametric problems.

    // Since we don't have a BlockOperator (yet? ;P), we assemble into one global matrix and use views for the local and
    // coupling matrices. We thus need a global DG DoF mapping.
    std::vector<size_t> index_offset(dd_grid.num_subdomains() + 1);
    size_t global_size = 0;
    for (size_t ss = 0; ss < dd_grid.num_subdomains(); ++ss) {
      index_offset[ss] = global_size;
      global_size += local_spaces[ss]->mapper().size();
    }
    index_offset[dd_grid.num_subdomains()] = global_size;

    // Assemble local sparsity patterns into a global one.
    XT::LA::SparsityPatternDefault global_pattern(global_size);
    const auto pattern_to_global = [&](const auto& local_pattern, const size_t ss, const size_t nn) {
      for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
        const auto global_ii = index_offset[ss] + local_ii;
        for (const auto& local_jj : local_pattern.inner(local_ii)) {
          const auto global_jj = index_offset[nn] + local_jj;
          global_pattern.insert(global_ii, global_jj);
        }
      }
    };
    for (auto&& macro_element : elements(macro_grid_view)) {
      const auto ss = dd_grid.subdomain(macro_element);
      const auto local_pattern = make_sparsity_pattern(*local_spaces[ss]);
      pattern_to_global(local_pattern, ss, ss);
      for (auto&& macro_intersection : intersections(macro_grid_view, macro_element)) {
        if (macro_intersection.neighbor()) {
          const auto& macro_neighbor = macro_intersection.outside();
          // Due to the nature of the coupling intersections, we don't have the hanging node problem. We can thus
          // treat each intersection only once.
          if (macro_grid_view.indexSet().index(macro_element) < macro_grid_view.indexSet().index(macro_neighbor)) {
            const size_t nn = dd_grid.subdomain(macro_neighbor);
            const auto coupling_pattern = make_coupling_sparsity_pattern(
                *local_spaces[ss],
                *local_spaces[nn],
                XT::Grid::make_coupling_grid_view(macro_element, macro_neighbor, dd_grid, macro_intersection));
            // The operator might be non-symmetric (for non-symmetric diffusion), but the sparsity pattern is.
            pattern_to_global(coupling_pattern, ss, nn);
            pattern_to_global(coupling_pattern, nn, ss);
          }
        }
      }
    }
    M global_matrix(global_size, global_size, global_pattern);
    V global_vector(global_size);

    logger.info() << "done" << std::endl;
    logger.info() << "Assembling subdomain contributions ... " << std::flush;

    // Assemble local discretization: any approximation of the original problem with zero-Neumann boundary values.
    for (auto&& macro_element : elements(macro_grid_view)) {
      const auto ss = dd_grid.subdomain(macro_element);
      const auto subdomain_grid = dd_grid.local_grid(ss).leaf_view();
      // LHS
      BilinearForm<GV_L> local_bilinear_form(subdomain_grid);
      local_bilinear_form += LocalElementIntegralBilinearForm<E_L>(LocalLaplaceIntegrand<E_L>(diffusion));
      if (!local_spaces[ss]->continuous(0)) {
        local_bilinear_form +=
            {LocalCouplingIntersectionIntegralBilinearForm<I_L>(
                 LocalLaplaceIPDGIntegrands::InnerCoupling<I_L>(symmetry_prefactor, diffusion, weight)
                 + LocalIPDGIntegrands::InnerPenalty<I_L>(penalty_parameter, weight)),
             XT::Grid::ApplyOn::InnerIntersectionsOnce<GV_L>()};
      }
      XT::LA::MatrixView<M> local_matrix(
          global_matrix, index_offset[ss], index_offset[ss + 1], index_offset[ss], index_offset[ss + 1]);
      auto subdomain_operator = local_bilinear_form.with(*local_spaces[ss], *local_spaces[ss], local_matrix);
      XT::LA::VectorView<V> local_vector(global_vector, index_offset[ss], index_offset[ss + 1]);
      // RHS
      auto subdomain_functional = make_vector_functional(*local_spaces[ss], local_vector);
      subdomain_functional.append(LocalElementIntegralFunctional<E_L>(LocalProductIntegrand<E_L>().with_ansatz(force)));
      // assembly
      auto subdomain_walker = XT::Grid::make_walker(subdomain_grid);
      subdomain_walker.append(subdomain_operator);
      subdomain_walker.append(subdomain_functional);
      subdomain_walker.walk(DXTC_CONFIG_GET("parallel", true));
    }

    logger.info() << "done" << std::endl;
    logger.info() << "Assembling coupling contributions ... " << std::flush;

    // Assemble coupling contributions.
    for (auto&& macro_element : elements(macro_grid_view)) {
      const auto ss = dd_grid.subdomain(macro_element);
      for (auto&& macro_intersection : intersections(macro_grid_view, macro_element)) {
        if (macro_intersection.neighbor()) {
          const auto& macro_neighbor = macro_intersection.outside();
          if (macro_grid_view.indexSet().index(macro_element) < macro_grid_view.indexSet().index(macro_neighbor)) {
            const size_t nn = dd_grid.subdomain(macro_neighbor);
            const auto coupling_grid =
                XT::Grid::make_coupling_grid_view(macro_element, macro_neighbor, dd_grid, macro_intersection);
            using GV_C = std::remove_const_t<decltype(coupling_grid)>;
            using I_C = XT::Grid::extract_intersection_t<GV_C>;
            // We need to fully specify the template arguments, to allow for GV_L != GV_C.
            BilinearForm<GV_C, 1, 1, 1, 1, double, GV_L, GV_L> coupling_bilinear_form(coupling_grid);
            // No need to restrict the intersection this local bilinear form is applied on, since the coupling grid only
            // contains the correct ones.
            coupling_bilinear_form += LocalCouplingIntersectionIntegralBilinearForm<I_C>(
                LocalLaplaceIPDGIntegrands::InnerCoupling<I_C>(symmetry_prefactor, diffusion, weight)
                + LocalIPDGIntegrands::InnerPenalty<I_C>(penalty_parameter, weight));
            // The above bilinear form has to be assembled into four matrices, corresponding to inside/outside and
            // test/ansatz contributions.
            // - we thus need need four matrices (note: rows = test 1st args, cols = ansatz 2nd args)
            XT::LA::MatrixView<M> coupling_matrix_ss_ss(
                global_matrix, index_offset[ss], index_offset[ss + 1], index_offset[ss], index_offset[ss + 1]);
            XT::LA::MatrixView<M> coupling_matrix_ss_nn(
                global_matrix, index_offset[ss], index_offset[ss + 1], index_offset[nn], index_offset[nn + 1]);
            XT::LA::MatrixView<M> coupling_matrix_nn_ss(
                global_matrix, index_offset[nn], index_offset[nn + 1], index_offset[ss], index_offset[ss + 1]);
            XT::LA::MatrixView<M> coupling_matrix_nn_nn(
                global_matrix, index_offset[nn], index_offset[nn + 1], index_offset[nn], index_offset[nn + 1]);
            // - and four operators (note: first ansatz space, then test space in make_matrix_operator!)
            auto coupling_operator_ss_ss =
                make_matrix_operator(coupling_grid, *local_spaces[ss], *local_spaces[ss], coupling_matrix_ss_ss);
            auto coupling_operator_ss_nn =
                make_matrix_operator(coupling_grid, *local_spaces[nn], *local_spaces[ss], coupling_matrix_ss_nn);
            auto coupling_operator_nn_ss =
                make_matrix_operator(coupling_grid, *local_spaces[ss], *local_spaces[nn], coupling_matrix_nn_ss);
            auto coupling_operator_nn_nn =
                make_matrix_operator(coupling_grid, *local_spaces[nn], *local_spaces[nn], coupling_matrix_nn_nn);
            // - and need to associate each operator with the respective contribution from the bilinear form
            coupling_operator_ss_ss.append(coupling_bilinear_form,
                                           {},
                                           {
                                               false, // element contribution
                                               true, //  coupling intersetcion contributions -  in/in
                                               false, //                                     -  in/out
                                               false, //                                     - out/in
                                               false, //                                     - out/out
                                               false //  non-coupling intersection contributions
                                           });
            coupling_operator_ss_nn.append(coupling_bilinear_form, {}, {false, false, true, false, false, false});
            coupling_operator_nn_ss.append(coupling_bilinear_form, {}, {false, false, false, true, false, false});
            coupling_operator_nn_nn.append(coupling_bilinear_form, {}, {false, false, false, false, true, false});
            // - and finally carry out the actual assembly
            auto coupling_walker = XT::Grid::make_walker(coupling_grid);
            coupling_walker.append(coupling_operator_ss_ss);
            coupling_walker.append(coupling_operator_ss_nn);
            coupling_walker.append(coupling_operator_nn_ss);
            coupling_walker.append(coupling_operator_nn_nn);
            coupling_walker.walk(DXTC_CONFIG_GET("parallel", true));
          }
        }
      }
    }

    // TODO: Assemble boundary contributions.

    // TODO: solve sglobal system and visualize solution
    //    auto solution = make_discrete_function<V>(space);
    //    XT::LA::make_solver(lhs_op.matrix()).apply(rhs_func.vector(), solution.dofs().vector());

    //    solution.visualize("solution");
  } catch (Exception& e) {
    std::cerr << "\nDUNE reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "\nstl reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occured!" << std::endl;
    return EXIT_FAILURE;
  } // try
  return EXIT_SUCCESS;
} // ... main(...)
