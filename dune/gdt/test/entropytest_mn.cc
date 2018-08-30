// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#include <iostream>

#include "config.h"

#if HAVE_TBB
#include <tbb/task_scheduler_init.h>
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/hyperbolic/mn-discretization.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/kinetictransport/testcases.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  // ***************** parse arguments and set up MPI and TBB
  size_t num_threads = 1;
  size_t threading_partition_factor = 1;
  size_t num_save_steps = 10;
  size_t num_output_steps = num_save_steps;
  size_t quad_refinements = 0;
  size_t quad_order = 31;
  std::string grid_size("100"), overlap_size("2");
  double t_end = 0;
  std::string filename;
  int parsed_ret = parse_momentmodel_arguments(argc,
                                               argv,
                                               num_threads,
                                               threading_partition_factor,
                                               num_save_steps,
                                               num_output_steps,
                                               quad_refinements,
                                               quad_order,
                                               grid_size,
                                               overlap_size,
                                               t_end,
                                               filename);
#if HAVE_TBB
  tbb::task_scheduler_init tbb_init(boost::numeric_cast<int>(num_threads));
#endif
  if (parsed_ret)
    return parsed_ret;

  static constexpr bool reconstruct = true;
  static constexpr size_t momentOrder = 0; // replacethisline

  using BasisfunctionType =
      // LegendreMomentBasis<double, double, momentOrder>; // 1d
      RealSphericalHarmonicsMomentBasis<double, double, momentOrder>; // 3d
  // HatFunctionMomentBasis<double, 1, double, momentOrder, 1, 1>; // 1d
  // HatFunctionMomentBasis<double, 3, double, momentOrder, 1, 3>; // 3d
  // PartialMomentBasis<double, 1, double, momentOrder, 1, 1>; // 1d
  // PartialMomentBasis<double, 3, double, momentOrder, 1, 3>; // 3d

  using TestCaseType =
      //    Hyperbolic::Problems::KineticTransport::PlaneSourceMnTestCase<Yasp1Grid, BasisfunctionType, reconstruct>;
      //    Hyperbolic::Problems::KineticTransport::SourceBeamMnTestCase<Yasp1Grid, BasisfunctionType, reconstruct>;
      //      Hyperbolic::Problems::KineticTransport::PointSourceMnTestCase<Yasp3Grid, BasisfunctionType, reconstruct>;
      Hyperbolic::Problems::KineticTransport::CheckerboardMnTestCase<Yasp3Grid, BasisfunctionType, reconstruct>;
  //      Hyperbolic::Problems::KineticTransport::ShadowMnTestCase<Yasp3Grid, BasisfunctionType, reconstruct>;

  HyperbolicMnDiscretization<TestCaseType> test;
  auto norms = test.run(
      num_save_steps, num_output_steps, quad_refinements, quad_order, grid_size, overlap_size, t_end, filename);
  std::cout << "l1norm = " << XT::Common::to_string(norms[0], 15) << std::endl;
  std::cout << "l2norm = " << XT::Common::to_string(norms[1], 15) << std::endl;
  std::cout << "linfnorm = " << XT::Common::to_string(norms[2], 15) << std::endl;
}
