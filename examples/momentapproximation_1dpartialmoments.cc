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

#include <dune/gdt/test/momentmodels/moment-approximation.hh>
#include <dune/gdt/test/momentmodels/basisfunctions.hh>


template <int momentOrder, Dune::GDT::EntropyType entropy>
struct moment_approximation_helper
{
  static void run(const std::string testcasename, const std::string filename)
  {
    using namespace Dune;
    using namespace Dune::GDT;

    using MomentBasisType = PartialMomentBasis<double, 1, double, momentOrder, 1, 1, 1, entropy>;

    using GridType = YASP_1D_EQUIDISTANT_OFFSET;
    using GridViewType = typename GridType::LeafGridView;
    using VectorType = typename Dune::XT::LA::Container<double, Dune::XT::LA::default_backend>::VectorType;
    using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType>;

    auto test = std::make_unique<MomentApproximation<MomentBasisType, DiscreteFunctionType>>();
    test->run(MomentBasisType::num_intervals, testcasename, filename);
    moment_approximation_helper<momentOrder - 2, entropy>::run(testcasename, filename);
  }
};

template <Dune::GDT::EntropyType entropy>
struct moment_approximation_helper<0, entropy>
{
  static void run(const std::string /*testcasename*/, const std::string /*filename*/) {}
};


int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::GDT;

  MPIHelper::instance(argc, argv);

  std::string testcasename = "Gauss1d";
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

  static constexpr int max_number_of_moments = 50;
  static_assert(!(max_number_of_moments % 2), "Maximal number of moments has to be even!");
  static constexpr EntropyType entropy = EntropyType::MaxwellBoltzmann;
  moment_approximation_helper<max_number_of_moments, entropy>::run(testcasename, testcasename);
}
