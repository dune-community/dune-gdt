// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include "config.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <utility>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/la/container/common/vector/dense.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/discretefunction/adaptation.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>

using namespace Dune;


// some global defines
using G = ONED_1D;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;
using V = XT::LA::CommonDenseVector<double>;


V compute_local_l2_norms(const XT::Functions::GridFunctionInterface<E>& func, const GV& grid_view)
{
  // the best index set for a grid view of arbitrary elements is a scalar FV space ...
  auto fv_space = GDT::make_finite_volume_space<1>(grid_view);
  // ... and the best way to associate data with grid elements is a corresponding discrete function
  auto df = GDT::make_discrete_function<V>(fv_space);
  // compute local L^2 products
  auto walker = XT::Grid::make_walker(grid_view);
  walker.append(/*prepare nothing*/ []() {},
                /*apply local*/
                [&](const auto& element) {
                  auto local_func = func.local_function();
                  local_func->bind(element);
                  // models \int_element 1*phi*psi dx for any phi/psi
                  const GDT::LocalElementIntegralBilinearForm<E> local_l2_product(
                      GDT::LocalElementProductIntegrand<E>(1.));
                  // evaluate local product with phi = psi = local_func
                  const auto element_l2_error2 = local_l2_product.apply2(*local_func, *local_func)[0][0];
                  // store in entry for this element (we keep them unsquared, makes for easier computation of total)
                  df.local_discrete_function(element)->dofs()[0] = std::sqrt(element_l2_error2);
                },
                /*finalize*/ []() {});
  walker.walk(/*thread_parallel=*/true);
  return df.dofs().vector();
} // ... compute_local_l2_norms(...)


void mark_elements(
    G& grid, const GV& grid_view, const V& local_indicators, const double& refine_factor, const double& coarsen_factor)
{
  // mark all elements with a contribution above coarsen_threashhold for coarsening
  // * therefore: as above, use a scalar FV space as index mapper
  auto fv_space = GDT::make_finite_volume_space<1>(grid_view);
  //              ... and a corresponding discrete function to associate data with grid elements
  //  auto local_indicators_function = GDT::make_discrete_function(fv_space, local_indicators);
  //  auto local_indicators_local_function = local_indicators_function.local_discrete_function();
  // prepare dörfler marking
  std::vector<std::pair<double, size_t>> accumulated_sorted_local_indicators(local_indicators.size());
  for (size_t ii = 0; ii < local_indicators.size(); ++ii)
    accumulated_sorted_local_indicators[ii] = {std::pow(local_indicators[ii], 2), ii};
  std::sort(accumulated_sorted_local_indicators.begin(),
            accumulated_sorted_local_indicators.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  for (size_t ii = 1; ii < accumulated_sorted_local_indicators.size(); ++ii)
    accumulated_sorted_local_indicators[ii].first += accumulated_sorted_local_indicators[ii - 1].first;
  // select smallest coarsen_factor contributions for coarsening
  std::set<size_t> elements_to_be_coarsened;
  const double total_indicators = std::pow(local_indicators.l2_norm(), 2);
  for (const auto& indicator : accumulated_sorted_local_indicators)
    if (indicator.first < (coarsen_factor * total_indicators))
      elements_to_be_coarsened.insert(indicator.second);
  // select largest refine_factor contributions
  std::set<size_t> elements_to_be_refined;
  for (const auto& indicator : accumulated_sorted_local_indicators)
    if (indicator.first > ((1 - refine_factor) * total_indicators))
      elements_to_be_refined.insert(indicator.second);
  // mark elements
  size_t corsend_elements = 0;
  size_t refined_elements = 0;
  for (auto&& element : elements(grid_view)) {
    const size_t index = fv_space.mapper().global_indices(element)[0];
    bool coarsened = false;
    if (std::find(elements_to_be_coarsened.begin(), elements_to_be_coarsened.end(), index)
        != elements_to_be_coarsened.end()) {
      grid.mark(/*coarsen*/ -1, element);
      coarsened = true;
    }
    if (std::find(elements_to_be_refined.begin(), elements_to_be_refined.end(), index)
        != elements_to_be_refined.end()) {
      grid.mark(/*refine and overwrite coarsening if present*/ 2, element);
      refined_elements += 1;
      coarsened = false;
    }
    if (coarsened)
      ++corsend_elements;
  }
  auto logger = XT::Common::TimedLogger().get("mark_elements");
  logger.info() << "marked " << corsend_elements << "/" << fv_space.mapper().size() << " elements for coarsening and "
                << refined_elements << "/" << fv_space.mapper().size() << " elements for refinement" << std::endl;
} // ... mark_elements(...)


int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // initial grid
    const double domain_length = 10;
    auto grid_ptr =
        XT::Grid::make_cube_grid<G>(/*lower_left=*/0., /*upper_right=*/domain_length, /*num_elements=*/1).grid_ptr();
    auto& grid = *grid_ptr;
    auto grid_view = grid.leafGridView();

    // funtion to interpolate
    const XT::Functions::GenericFunction<d> cosh(/*approx_pol_order=*/10,
                                                 [](const auto& x, const auto& /*param*/) { return std::cosh(x); });
    const auto reference_norm = compute_local_l2_norms(cosh.as_grid_function(grid_view), grid_view).l2_norm();

    // fake solution to demonstrate restriction/prolongation
    auto fv_space = GDT::make_finite_volume_space<1>(grid_view);
    auto current_solution = GDT::make_discrete_function<V>(fv_space);
    GDT::interpolate(cosh, current_solution);

    // the main adaptation loop
    auto helper = make_adaptation_helper(grid, fv_space, current_solution);
    const double tolerance = DXTC_CONFIG_GET("tolerance", 1e-3);
    size_t counter = 0;
    while (true) {
      logger.info() << "step " << counter << ", space has " << fv_space.mapper().size() << " DoFs" << std::endl;
      current_solution.visualize("interpolated_function_" + XT::Common::to_string(counter));

      // compute local L^2 errors
      const auto local_errors = compute_local_l2_norms(current_solution - cosh.as_grid_function(grid_view), grid_view);
      logger.info() << "  relative interpolation error: " << local_errors.l2_norm() / reference_norm << std::endl;
      if (local_errors.l2_norm() / reference_norm < tolerance) {
        logger.info() << "target accuracy reached, terminating!" << std::endl;
        break;
      }

      // use these as indicators for grid refinement
      mark_elements(grid,
                    grid_view,
                    local_errors,
                    DXTC_CONFIG_GET("mark.refine_factor", 0.25),
                    DXTC_CONFIG_GET("mark.coarsen_factor", 0.01));

      // adapt the grid (restricts/prolongs the solution)
      helper.pre_adapt();
      helper.adapt();
      helper.post_adapt();

      // now that the grid is adapted, interpolate the function on the new grid
      GDT::interpolate(cosh, current_solution);

      ++counter;
    }

    double min_h = std::numeric_limits<double>::max();
    for (auto&& element : elements(grid_view))
      min_h = std::min(min_h, XT::Grid::diameter(element));

    logger.info() << "\nA uniformly refinement grid of similar accuracy would have roughly required "
                  << int(domain_length / min_h) << " DoFs" << std::endl;

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
