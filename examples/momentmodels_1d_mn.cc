// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2019)

#include <iostream>

#include "config.h"

#if HAVE_TBB
#  include <tbb/task_scheduler_init.h>
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/mn-discretization.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  /****************** parse arguments and set up MPI and TBB *************************/

  size_t num_threads, threading_partition_factor, num_save_steps, num_output_steps, quad_order, quad_refinements,
      overlap_size;
  double t_end = 0;
  std::string filename, grid_size;
  parse_momentmodel_arguments(argc,
                              argv,
                              num_threads,
                              threading_partition_factor,
                              num_save_steps,
                              num_output_steps,
                              quad_order,
                              quad_refinements,
                              grid_size,
                              overlap_size,
                              t_end,
                              filename);
#if HAVE_TBB
  tbb::task_scheduler_init tbb_init(boost::numeric_cast<int>(num_threads));
#endif

  /***************** choose test case and approximation parameters  ***********************/

  // if reconstruct = false, a first order scheme without reconstruction is used
  static constexpr bool reconstruct = true;
  // moment_order is the polynomial order of the Legendre basis, has to be at least 1. The number of moments is
  // moment_order + 1.
  static constexpr size_t moment_order = 9;

  using MomentBasis = LegendreMomentBasis<double, double, moment_order>;

  // choose test case here
  using TestCaseType = PlaneSourceMnTestCase<YASP_1D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;
  // using TestCaseType = SourceBeamMnTestCase<YASP_1D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;

  /************************************ run **********************************************/

  HyperbolicMnDiscretization<TestCaseType> test;
  const auto norms_and_rank = test.run(
      num_save_steps, num_output_steps, quad_order, quad_refinements, grid_size, overlap_size, t_end, filename);
  const auto& norms = norms_and_rank.first;
  if (norms_and_rank.second == 0) {
    std::cout << "l1norm = " << XT::Common::to_string(norms[0], 15) << std::endl;
    std::cout << "l2norm = " << XT::Common::to_string(norms[1], 15) << std::endl;
    std::cout << "linfnorm = " << XT::Common::to_string(norms[2], 15) << std::endl;
  }
}
