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
  // num_moments is the number of moments, has to be at least 2 and can only take even values.
  static constexpr size_t num_moments = 10;

  using MomentBasis = PartialMomentBasis<double, 1, double, num_moments, 1, 1>;

  // choose test case here
  using TestCaseType = PlaneSourcePnTestCase<YASP_1D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;
  // using TestCaseType = SourceBeamPnTestCase<YASP_1D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;

  /************************************ run **********************************************/

  HyperbolicPnDiscretization<TestCaseType> test;
  const auto norms = test.run(
      num_save_steps, num_output_steps, quad_order, quad_refinements, grid_size, overlap_size, t_end, filename);
  std::cout << "l1norm = " << XT::Common::to_string(norms[0], 15) << std::endl;
  std::cout << "l2norm = " << XT::Common::to_string(norms[1], 15) << std::endl;
  std::cout << "linfnorm = " << XT::Common::to_string(norms[2], 15) << std::endl;
}
