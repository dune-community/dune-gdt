// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#include <iostream>

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/common/parallel/threadmanager.hh>

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/momentmodels/moment-approximation.hh>
#include <dune/gdt/momentmodels/basisfunctions.hh>


template <int refinement, Dune::GDT::EntropyType entropy>
struct moment_approximation_helper
{
  static void run(const int quadrature_refinements, const std::string testcasename, const std::string filename)
  {
    using namespace Dune;
    using namespace Dune::GDT;

    using BasisfunctionType = HatFunctionMomentBasis<double, 3, double, refinement, 1, 3, entropy>;

    using GridType = YASP_3D_EQUIDISTANT_OFFSET;
    using GridViewType = typename GridType::LeafGridView;
    using VectorType = typename Dune::XT::LA::Container<double, Dune::XT::LA::default_backend>::VectorType;
    using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType>;
    auto test = std::make_unique<MomentApproximation<BasisfunctionType, DiscreteFunctionType>>();
    test->run(quadrature_refinements, testcasename, filename);
    moment_approximation_helper<refinement - 1, entropy>::run(quadrature_refinements, testcasename, filename);
  }
};

template <Dune::GDT::EntropyType entropy>
struct moment_approximation_helper<-1, entropy>
{
  static void
  run(const int /*quadrature_refinements*/, const std::string /*testcasename*/, const std::string /*filename*/)
  {}
};


int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  MPIHelper::instance(argc, argv);

  std::string testcasename = "GaussOnSphere";
  if (argc >= 2)
    testcasename = argv[1];
  if (argc == 3) {
    DXTC_CONFIG["threading.max_count"] = argv[2];
    XT::Common::threadManager().set_max_threads(XT::Common::from_string<size_t>(argv[2]));
  } else if (argc > 3) {
    std::cerr << "Too many command line arguments, please provide a testcase name and the number of threads only!"
              << std::endl;
    return 1;
  }

  static constexpr EntropyType entropy = EntropyType::MaxwellBoltzmann;
  // static constexpr EntropyType entropy = EntropyType::BoseEinstein;
  static const int max_refinements = 4;
  const int quadrature_refinements = 6;
  if (quadrature_refinements < max_refinements)
    DUNE_THROW(Dune::InvalidStateException,
               "The quadrature has to use at least as many spherical triangles as the highest-order model!");
  moment_approximation_helper<max_refinements, entropy>::run(quadrature_refinements, testcasename, testcasename);
}
