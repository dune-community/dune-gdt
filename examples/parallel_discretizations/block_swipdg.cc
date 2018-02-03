// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk (2017 - 2018)

#define ALUGRIDDEBUG 1

#include <config.h>

#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif
#include <dune/alugrid/common/structuredgridfactory.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/parallel/threadmanager.hh>


#if HAVE_TBB
#include <tbb/task_scheduler_init.h>
#endif

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/common/misc.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/common/timings.hh>

#include <dune/gdt/projections.hh>

#include <dune/gdt/test/linearelliptic/discretizers/block-ipdg.hh>
#include <dune/gdt/test/linearelliptic/discretizers/cg.hh>
#include <dune/gdt/test/linearelliptic/problems/ESV2007.hh>
#include <dune/gdt/test/linearelliptic/problems/ER2007.hh>
#include <dune/gdt/operators/laplace.hh>

using namespace Dune;

static_assert(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS, "grid extensions are mandatory");

//  using GridType = Dune::ALUGrid<3, 3, simplex, nonconforming, ALUGridMPIComm>;
// using GridType = Dune::ALUGrid<2, 2, simplex, nonconforming, ALUGridMPIComm>;
// using GridType = Dune::ALUGrid<2, 2, cube, nonconforming, ALUGridMPIComm>;
using GridType = Yasp2Grid;
using GridProvider = XT::Grid::GridProvider<GridType, XT::Grid::DD::SubdomainGrid<GridType>>;
using Errors = std::pair<double, double>;
constexpr const size_t pol_order = 1;
constexpr const auto la_backend = XT::LA::Backends::istl_sparse;

template <GDT::Backends space_backend>
struct Traits
{
  using CGFactory =
      GDT::LinearElliptic::CGDiscretizer<GridType, XT::Grid::Layers::leaf, space_backend, la_backend, pol_order>;
  using CGDiscretization = typename CGFactory::DiscretizationType;
  using IpdgFactory =
      GDT::LinearElliptic::IpdgDiscretizer<GridType, XT::Grid::Layers::leaf, space_backend, la_backend, pol_order>;
  using IpdgDiscretization = typename IpdgFactory::DiscretizationType;

  static const auto local_method = GDT::LocalEllipticIpdgIntegrands::Method::swipdg;
  static const auto coupling_method = GDT::LocalEllipticIpdgIntegrands::Method::swipdg;
  using BlockIpdgFactory = GDT::LinearElliptic::
      BlockIpdgDiscretizer<GridType, space_backend, la_backend, pol_order, double, 1, local_method, coupling_method>;
  using BlockIpdgDiscretization = typename BlockIpdgFactory::DiscretizationType;
  using Vector = typename BlockIpdgDiscretization::VectorType;
  static_assert(std::is_same<Vector, typename IpdgDiscretization::VectorType>::value, "Vector mismatch");
  using IpdgFunction = GDT::DiscreteFunction<typename IpdgDiscretization::AnsatzSpaceType, Vector>;
  using CGFunction = GDT::DiscreteFunction<typename CGDiscretization::AnsatzSpaceType, Vector>;
  using BlockIpdgFunction = GDT::DiscreteFunction<typename BlockIpdgDiscretization::AnsatzSpaceType, Vector>;
  using ConstBlockIpdgFunction = GDT::ConstDiscreteFunction<typename BlockIpdgDiscretization::AnsatzSpaceType, Vector>;
};

std::string path(std::string filename)
{
  boost::filesystem::path ret(DXTC_CONFIG_GET("global.datadir", "data"));
  return (ret / filename).string();
}

template <GDT::Backends spc, class Space, class TestCase>
Errors exact_error(const TestCase& test_case,
                   const GDT::ConstDiscreteFunction<Space, typename Traits<spc>::Vector>& discrete_solution)
{
  XT::Common::ScopedTiming timing(std::string("analytical_error.") + discrete_solution.name());
  auto layer = discrete_solution.space().grid_layer();
  const auto& exact = test_case.exact_solution();
  auto diff = discrete_solution - exact;
  if (DXTC_CONFIG_GET("global.visualize", false)) {
    diff.visualize(layer, path(discrete_solution.name() + "_diff"), "solution");
  }
  const auto l2_error = Dune::GDT::make_l2_operator(layer, 2)->induced_norm(diff);
  const auto h1s_error = Dune::GDT::make_laplace_operator(layer, 2)->induced_norm(diff);
  return {l2_error, h1s_error};
}

template <GDT::Backends spc, class TestCase>
std::pair<typename Traits<spc>::Vector, typename Traits<spc>::CGDiscretization::AnsatzSpaceType>
cg(TestCase& test_case, GridProvider& grid_provider)
{
  XT::Common::ScopedTiming timer_all("cg.all");
  XT::Common::timings().start("cg.discretize");
  auto cg_disc = Traits<spc>::CGFactory::discretize(grid_provider, test_case.problem());
  auto cg_solution_vector = cg_disc.create_vector();
  XT::Common::timings().stop("cg.discretize");
  const auto solver_options = DXTC_CONFIG.sub("solver");
  XT::Common::timings().start("cg.solve");
  try {
    XT::Common::ScopedTiming solve("cg.solve");
    cg_disc.solve(cg_solution_vector, solver_options);
  } catch (XT::LA::Exceptions::linear_solver_failed_bc_it_did_not_converge) {
    cg_solution_vector += std::nan("failed");
  }
  if (DXTC_CONFIG_GET("global.visualize", false)) {
    cg_disc.visualize(cg_solution_vector, path("cg_solution"), "solution");
  }
  return {cg_solution_vector, cg_disc.ansatz_space()};
}

template <GDT::Backends spc, class BlockTestCase>
std::pair<typename Traits<spc>::Vector, typename Traits<spc>::BlockIpdgDiscretization::AnsatzSpaceType>
block_ipdg(BlockTestCase& test_case, GridProvider& grid_provider)
{
  XT::Common::ScopedTiming timer_all("block_ipdg.all");
  XT::Common::timings().start("block_ipdg.discretize");
  auto block_ipdg_disc = Traits<spc>::BlockIpdgFactory::discretize(grid_provider, test_case.problem(), /*not used*/ 0);
  auto block_ipdg_solution_vector = block_ipdg_disc.create_vector();
  XT::Common::timings().stop("block_ipdg.discretize");
  const auto solver_options = DXTC_CONFIG.sub("solver");
  XT::Common::timings().start("block_ipdg.solve");
  try {
    XT::Common::ScopedTiming solve("block_ipdg.solve");
    block_ipdg_disc.solve(block_ipdg_solution_vector, solver_options);
  } catch (XT::LA::Exceptions::linear_solver_failed_bc_it_did_not_converge) {
    block_ipdg_solution_vector += std::nan("failed");
  }
  if (DXTC_CONFIG_GET("global.visualize", false)) {
    block_ipdg_disc.visualize(block_ipdg_solution_vector, path("block_ipdg_solution"), "solution");
  }
  return {block_ipdg_solution_vector, block_ipdg_disc.ansatz_space()};
}

template <GDT::Backends spc, class TestCase>
std::pair<typename Traits<spc>::IpdgDiscretization, typename Traits<spc>::Vector>
plain_ipdg(TestCase& test_case, GridProvider& grid_provider)
{
  XT::Common::ScopedTiming timer_all("ipdg.all");
  XT::Common::timings().start("ipdg.discretize");
  auto ipdg_disc = Traits<spc>::IpdgFactory::discretize(grid_provider, test_case.problem());
  auto ipdg_solution_vector = ipdg_disc.create_vector();
  const auto solver_options = DXTC_CONFIG.sub("solver");
  XT::Common::timings().stop("ipdg.discretize");
  try {
    XT::Common::ScopedTiming solve("ipdg.solve");
    ipdg_disc.solve(ipdg_solution_vector, solver_options);
  } catch (XT::LA::Exceptions::linear_solver_failed_bc_it_did_not_converge) {
    ipdg_solution_vector += std::nan("failed");
  }
  if (DXTC_CONFIG_GET("global.visualize", false)) {
    ipdg_disc.visualize(ipdg_solution_vector, path("ipdg_solution"), "solution");
  }
  return {std::move(ipdg_disc), std::move(ipdg_solution_vector)};
}

template <GDT::Backends spc, class TestCase>
void error_output(const std::map<std::string, Errors>& error_map, const TestCase& test_case)
{
  auto& out = DXTC_LOG_INFO_0;
  const auto spc_str = GDT::backend_names[spc];
  const auto prob_str = XT::Common::get_template_basename(test_case.problem());
  out << "\n"
      << prob_str << " | " << spc_str << "\n"
      << "----------------------------------\n";
  boost::format row_fmt("%8s %|8t| %15.9f %|8t|%15.9f\n");
  boost::format head_fmt("%8s %|8t| %-15s %|8t|%-15s\n");
  // header
  out << head_fmt % "" % "L2" % "H1s";
  for (auto&& row : error_map) {
    const auto disc = row.first;
    const auto errors = row.second;
    out << row_fmt % disc % errors.first % errors.second;
  }

  out.flush();
}

template <GDT::Backends spc, class TestCase, class BlockTestCase>
void single_run()
{
  XT::Common::ScopedTiming timer_all("single_run");

  std::map<std::string, Errors> error_map;
  XT::Common::Configuration cfg(TestCase::grid_cfg(), std::string("grid"));
  XT::Common::Configuration block_cfg(BlockTestCase::grid_cfg(), std::string("grid"));

  //! get and overwrite grid values from global config
  block_cfg.add(DXTC_CONFIG, "", true);
  cfg.add(DXTC_CONFIG, "", true);

  // Testcases are already grid providers and grid privder can return grid provider
  TestCase test_case{cfg};
  BlockTestCase block_test_case{block_cfg};
  auto& grid_provider = block_test_case.level_provider(0);

  auto ipdg_disc_and_solution = plain_ipdg<spc>(test_case, grid_provider);
  const auto& space = ipdg_disc_and_solution.first.ansatz_space();
  const typename Traits<spc>::IpdgFunction ipdg_solution(space, ipdg_disc_and_solution.second, "ipdg");

  if (DXTC_CONFIG_GET("global.ipdg", true)) {
    error_map["ipdg"] = exact_error<spc>(test_case, ipdg_solution);
  }

  if (DXTC_CONFIG_GET("global.cg", true)) {
    auto cg_disc_and_solution = cg<spc>(test_case, grid_provider);
    typename Traits<spc>::CGFunction cg_func(cg_disc_and_solution.second, cg_disc_and_solution.first, "cg");
    error_map["cg"] = exact_error<spc>(test_case, cg_func);
  }

  if (DXTC_CONFIG_GET("global.block", false)) {
    const auto block_solution_vector_and_space = block_ipdg<spc>(block_test_case, grid_provider);
    typename Traits<spc>::ConstBlockIpdgFunction block_ipdg_solution(
        block_solution_vector_and_space.second, block_solution_vector_and_space.first, "block_ipdg");
    error_map["block"] = exact_error<spc>(test_case, block_ipdg_solution);
  }
  error_output<spc>(error_map, test_case);
}

int main(int argc, char** argv)
{
  using namespace Dune::XT::Common;
  try {

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    DXTC_CONFIG.read_command_line(argc, argv);

    test_create_directory(DXTC_CONFIG_GET("global.datadir", "data/"));

    // LOG_NONE = 1, LOG_ERROR = 2, LOG_INFO = 4,LOG_DEBUG = 8,LOG_CONSOLE =
    // 16,LOG_FILE = 32
    // --> LOG_ERROR | LOG_INFO | LOG_DEBUG | LOG_CONSOLE | LOG_FILE = 62
    Logger().create(DXTC_CONFIG_GET("logging.level", 62),
                    DXTC_CONFIG_GET("logging.file", std::string(argv[0]) + ".log"),
                    DXTC_CONFIG_GET("global.datadir", "data"),
                    DXTC_CONFIG_GET("logging.dir", "log" /*path below datadir*/));
    const int max_level = DXTC_CONFIG_GET("logging.max_nesting", -1);
    //    TimedLogger().create(max_level /*info*/, max_level /*debug*/);
    DXTC_TIMINGS.set_outputdir(DXTC_CONFIG_GET("global.datadir", "data"));

    //    DXTC_LOG_INFO_0 << DXTC_CONFIG.report_string() << std::endl;

    threadManager().set_max_threads(1);
    {
      OutputScopedTiming outs("all", DXTC_LOG_INFO_0);
      DXTC_CONFIG.set("grids.total_macro_cells", 256);

      switch (DXTC_CONFIG_GET("global.problem", 0)) {
        case 0:
          single_run<GDT::Backends::gdt,
                     GDT::LinearElliptic::ER2007TestCase<GridType>,
                     GDT::LinearElliptic::ER2007DdSubdomainsTestCase<GridType>>();
          break;
        case 1:
          single_run<GDT::Backends::gdt,
                     GDT::LinearElliptic::ESV2007TestCase<GridType>,
                     GDT::LinearElliptic::ESV2007DdSubdomainsTestCase<GridType>>();
          break;
        case 2:
          single_run<GDT::Backends::gdt,
                     GDT::LinearElliptic::ER2007TestCase<GridType>,
                     GDT::LinearElliptic::ER2007DdSubdomainsTestCase<GridType>>();
          single_run<GDT::Backends::gdt,
                     GDT::LinearElliptic::ESV2007TestCase<GridType>,
                     GDT::LinearElliptic::ESV2007DdSubdomainsTestCase<GridType>>();
          break;
      }
    }
    DXTC_TIMINGS.output_per_rank("profiler");
    mem_usage();
    //    dump_environment();

  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  }
}
