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

#include <dune/gdt/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/mn-discretization.hh>
#include <dune/gdt/test/momentmodels/kinetictransport/testcases.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  /****************** parse arguments and set up MPI and TBB *************************/

  size_t num_threads, threading_partition_factor, num_save_steps, num_output_steps, quad_order, quad_refinements,
      grid_size, overlap_size;
  double t_end = 0;
  std::string filename;
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
  // num_refinements is the number of times the initial triangulation of the sphere is refined. We start with the
  // octaeder with corners [+-1, 0, 0], [0, +-1, 0], [0, 0, +-1] projected to the sphere, thus there are 8 initial
  // triangles. At each refinement, each of the triangles is divided in 4 smaller triangles. For example, after one
  // refinement, there are 32 triangles. There are 4 partial moment basis functions per triangle and thus
  // 8*4^{num_refinements+1} moments.
  static constexpr size_t num_refinements = 0;

  using MomentBasis = PartialMomentBasis<double, 3, double, num_refinements, 1, 3>;

  // choose test case here
  using TestCaseType = PointSourceMnTestCase<YASP_3D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;
  // using TestCaseType = CheckerboardMnTestCase<YASP_3D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;
  // using TestCaseType = ShadowMnTestCase<YASP_3D_EQUIDISTANT_OFFSET, MomentBasis, reconstruct>;

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
