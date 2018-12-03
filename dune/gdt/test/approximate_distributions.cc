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

#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/hyperbolic/moment-approximation.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions.hh>

template <int momentOrder>
struct moment_approximation_helper
{
  static void run(const size_t quad_refinements, const size_t quad_order, const std::string filename)
  {
    using namespace Dune;
    using namespace Dune::GDT;

    using BasisfunctionType =
        //         LegendreMomentBasis<double, double, momentOrder>; // 1d
        //        RealSphericalHarmonicsMomentBasis<double, double, momentOrder>; // 3d
        //    HatFunctionMomentBasis<double, 1, double, momentOrder, 1, 1>; // 1d
        HatFunctionMomentBasis<double, 3, double, momentOrder, 1, 3>; // 3d
    //     PartialMomentBasis<double, 1, double, momentOrder, 1, 1>; // 1d
    //     PartialMomentBasis<double, 3, double, momentOrder, 1, 3>; // 3d


    using GridType = Yasp3Grid;
    using GridLayerType = typename GridType::LeafGridView;
    using SpaceType = FvSpace<GridLayerType, double, BasisfunctionType::dimRange, 1>;
    using VectorType = typename Dune::XT::LA::Container<double, Dune::XT::LA::default_backend>::VectorType;
    using DiscreteFunctionType = DiscreteFunction<SpaceType, VectorType>;

    auto test = std::make_unique<MomentApproximation<BasisfunctionType, DiscreteFunctionType>>();
    test->run(quad_refinements - momentOrder, quad_order, filename);
    //    test.run(quad_refinements, quad_order, filename);
    moment_approximation_helper<momentOrder - 1>::run(quad_refinements, quad_order, filename);
  }
};

template <>
// struct moment_approximation_helper<0>
struct moment_approximation_helper<-1>
{
  static void run(const size_t /*quad_refinements*/, const size_t /*quad_order*/, const std::string /*filename*/) {}
};


int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;
  MPIHelper::instance(argc, argv);
#if HAVE_TBB
  DXTC_CONFIG.set("threading.partition_factor", 10, true);
  XT::Common::threadManager().set_max_threads(32);
#endif

  // ***************** parse arguments and set up MPI and TBB
  size_t quad_refinements = 0;
  size_t quad_order = 197;
  std::string filename;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--filename") {
      if (i + 1 < argc) {
        filename = std::string(argv[++i]);
      } else {
        std::cerr << "--filename option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--quad_refinements") {
      if (i + 1 < argc) {
        quad_refinements = XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "--quad_refinements option requires one argument." << std::endl;
        return 1;
      }
    } else if (std::string(argv[i]) == "--quad_order") {
      if (i + 1 < argc) {
        quad_order = XT::Common::from_string<size_t>(argv[++i]);
      } else {
        std::cerr << "--quad_order option requires one argument." << std::endl;
        return 1;
      }
    } else {
      std::cerr << "Unknown option " << std::string(argv[i]) << std::endl;
      return 1;
    }
  }

  // moment_approximation_helper<10>::run(quad_refinements, quad_order, filename);
  moment_approximation_helper<2>::run(quad_refinements, quad_order, filename);
}
