// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <examples/coordinate-transformed-mn.hh>

#include <chrono>
#include <random>
#include <thread>

using namespace Dune;
using namespace Dune::GDT;

template <class ModelSolver>
void check(ModelSolver& model_solver, const size_t num_output_dofs)
{
  using VectorType = typename ModelSolver::VectorType;
  // random number generator
  std::random_device dev;
  std::mt19937 rng(dev());
  // initial_values
  std::uniform_real_distribution<double> double_distrib(-1., 1.);
  const auto initial_values = model_solver.get_initial_values();
  const size_t num_dofs = initial_values.size();
  // const size_t num_output_dofs = num_dofs;
  std::vector<size_t> output_dofs(num_output_dofs);
  std::uniform_int_distribution<typename std::mt19937::result_type> distrib(0, num_dofs - 1);
  auto source_vec = initial_values;
  // ensure source_vec entries are non-zero to avoid masking errors
  for (auto& entry : source_vec)
    entry += double_distrib(rng);
  std::chrono::duration<double> restricted_prep_time(0.);
  std::chrono::duration<double> restricted_apply_time(0.);
  std::chrono::duration<double> prep_time(0.);
  std::chrono::duration<double> apply_time(0.);
  for (size_t run = 0; run < 10; ++run) {
    std::cout << "1d run " << run << std::endl;
    for (size_t ii = 0; ii < num_output_dofs; ++ii)
      output_dofs[ii] = distrib(rng);
    // std::sort(output_dofs.begin(), output_dofs.end());
    std::cout << "Restricted start" << std::endl;
    auto begin = std::chrono::steady_clock::now();
    model_solver.prepare_restricted_operator(output_dofs);
    restricted_prep_time += std::chrono::steady_clock::now() - begin;
    const auto& source_dofs = model_solver.restricted_op_input_dofs();
    const size_t num_source_dofs = source_dofs.size();
    VectorType restricted_source(num_source_dofs, 0.);
    for (size_t ii = 0; ii < num_source_dofs; ++ii)
      restricted_source[ii] = source_vec[source_dofs[ii]];
    begin = std::chrono::steady_clock::now();
    auto restricted_result = model_solver.apply_restricted_operator(restricted_source);
    restricted_apply_time += std::chrono::steady_clock::now() - begin;
    begin = std::chrono::steady_clock::now();
    auto result = model_solver.apply_operator(source_vec);
    apply_time += std::chrono::steady_clock::now() - begin;
    const double apply_tol = 1e-14;
    for (size_t ii = 0; ii < num_output_dofs; ++ii) {
      if (XT::Common::FloatCmp::ne(restricted_result[ii], result[output_dofs[ii]], apply_tol, apply_tol))
        std::cout << "Failed apply restricted: " << ii << ", " << output_dofs[ii] << ", " << result[output_dofs[ii]]
                  << ", " << restricted_result[ii]
                  << ", diff: " << XT::Common::to_string(std::abs(result[output_dofs[ii]] - restricted_result[ii]), 15)
                  << std::endl;
    } // ii
  } // runs
  std::cout << "prep: " << prep_time.count() << "  vs. " << restricted_prep_time.count() << std::endl;
  std::cout << "apply: " << apply_time.count() << "  vs. " << restricted_apply_time.count() << std::endl;
} // void check(...)

int main(int argc, char* argv[])
{
  using namespace Dune;
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
#if HAVE_TBB
    DXTC_CONFIG.set("threading.partition_factor", 1, true);
    XT::Common::threadManager().set_max_threads(1);
#endif

    using GridType1d = YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>;
    using GridType3d = YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>;
    using GV1d = typename GridType1d::LeafGridView;
    using GV3d = typename GridType3d::LeafGridView;
    using HF50Basis1d = HatFunctionMomentBasis<double, 1, double, 50>;
    // using PM50Basis1d = PartialMomentBasis<double, 1, double, 50>;
    // using M50Basis1d = LegendreMomentBasis<double, double, 50>;
    // using HF66Basis3d = HatFunctionMomentBasis<double, 3, double, /*refinements = */ 2>;
    // using PM128Basis3d = PartialMomentBasis<double, 3, double, /*refinements = */ 1>;
    using M3Basis3d = RealSphericalHarmonicsMomentBasis<double, double, /*order = */ 3, /*fluxdim = */ 3>;
    using HFM50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, HF50Basis1d>>;
    // using PMM50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, PM50Basis1d>>;
    // using M50SourceBeamSolver = CoordinateTransformedMnSolver<SourceBeamMn<GV1d, M50Basis1d>>;
    // using HFM66CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, HF66Basis3d>>;
    // using PMM128CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, PM128Basis3d>>;
    using M3CheckerboardSolver3d = CoordinateTransformedMnSolver<CheckerboardMn<GV3d, M3Basis3d>>;
    using ModelSolver1d = HFM50SourceBeamSolver;
    using ModelSolver3d = M3CheckerboardSolver3d;

    ModelSolver1d model_solver_1d("1d", /*num_save_steps*/ 1, /*grid_size*/ 1000);
    const size_t num_output_dofs_1d = 50;
    check(model_solver_1d, num_output_dofs_1d);

    ModelSolver3d model_solver_3d("3d", 1, 21);
    const size_t num_output_dofs_3d = 50;
    check(model_solver_3d, num_output_dofs_3d);

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
