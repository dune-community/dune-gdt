// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_AO2013_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_AO2013_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/functions/checkerboard.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/stationary-testcase.hh>
#include <dune/gdt/test/grids.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template <class E, class D, int d, class R, int r = 1>
class AO2013Problem : public ProblemBase<E, D, d, R, r>
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class AO2013Problem<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
    : public ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ScalarConstantFunctionType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> MatrixConstantFunctionType;
  typedef Dune::XT::Functions::CheckerboardFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
      CheckerboardFunctionType;

public:
  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration cfg = XT::Grid::cube_gridprovider_default_config();
    cfg["lower_left"] = "[0 0]";
    cfg["upper_right"] = "[1 1]";
    return cfg;
  }

  static XT::Common::Configuration default_boundary_info_cfg()
  {
    return XT::Grid::alldirichlet_boundaryinfo_default_config();
  }

  AO2013Problem(const XT::Common::Configuration& grd_cfg = default_grid_cfg(),
                const XT::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(
          new CheckerboardFunctionType({0.0, 0.0},
                                       {1.0, 1.0},
                                       {6, 6},
                                       {1.0, 1.0,  1.0, 0.1, 0.1, 0.1, 1.0, 0.01, 1.0, 0.1, 0.1, 0.1,
                                        1.0, 1.0,  1.0, 0.1, 0.1, 0.1, 1.0, 1.0,  1.0, 0.1, 0.1, 0.1,
                                        1.0, 0.01, 1.0, 0.1, 0.1, 0.1, 1.0, 1.0,  1.0, 0.1, 0.1, 0.1},
                                       "diffusion_factor"),
          new MatrixConstantFunctionType(XT::Functions::internal::unit_matrix<RangeFieldImp, 2>(), "diffusion_tensor"),
          new ScalarConstantFunctionType(1, "force"),
          new ScalarConstantFunctionType(0, "dirichlet"),
          new ScalarConstantFunctionType(0, "neumann"),
          grd_cfg,
          bnd_cfg)
  {
  }
}; // class AO2013Problem< ..., 1 >


template <class G, class R = double, int r = 1>
class AO2013TestCase
    : public Test::StationaryTestCase<G,
                                      LinearElliptic::AO2013Problem<typename G::template Codim<0>::Entity,
                                                                    typename G::ctype,
                                                                    G::dimension,
                                                                    R,
                                                                    r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  typedef LinearElliptic::AO2013Problem<E, D, d, R, r> ProblemType;

private:
  typedef Test::StationaryTestCase<G, ProblemType> BaseType;

#if DXT_DISABLE_LARGE_TESTS

  template <class T, bool anything = true>
  struct Helper
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[4 4]";
      return cfg;
    }
  };

#else // DXT_DISABLE_LARGE_TESTS

  template <class T, bool anything = true>
  struct Helper
  {
    static_assert(AlwaysFalse<T>::value, "Please add a configuration for this grid type!");
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<Yasp2Grid, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[6 6]";
      return cfg;
    }
  };

#if HAVE_DUNE_ALUGRID
  template <bool anything>
  struct Helper<AluConform2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[6 6]";
      cfg["num_refinements"] = "1";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluCube2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[6 6]";
      cfg["num_refinements"] = "1";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluSimplex2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[6 6]";
      return cfg;
    }
  };

#endif // HAVE_DUNE_ALUGRID
#endif // DXT_DISABLE_LARGE_TESTS

  static XT::Common::Configuration grid_cfg()
  {
    auto cfg = ProblemType::default_grid_cfg();
    cfg = Helper<typename std::decay<G>::type>::value(cfg);
    return cfg;
  }

public:
  using typename BaseType::GridType;

  AO2013TestCase(const size_t num_refs =
#if DXT_DISABLE_LARGE_TESTS
                     1
#else
                     3
#endif
                 )
    : BaseType(grid_cfg(), num_refs)
    , problem_()
  {
  }

  AO2013TestCase(XT::Common::Configuration cfg)
    : BaseType(cfg.sub("grid", false, grid_cfg()), cfg.get("grid.num_refinements", 1))
    , problem_()
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual void print_header(std::ostream& out = DXTC_LOG_INFO_0) const override final
  {
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: AO2013, local thermal block problem                     ||\n"
        << "||  (see http://wwwmath.uni-muenster.de/num/publications/2013/AO13/)  ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                                          ||\n"
        << "||  diffusion:  see page 3 (mu_test)                                  ||\n"
        << "||  force     = 1                                                     ||\n"
        << "||  dirichlet = 0                                                     ||\n"
        << "||  reference solution: discrete solution on finest grid              ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class AO2013TestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_AO2013_HH
