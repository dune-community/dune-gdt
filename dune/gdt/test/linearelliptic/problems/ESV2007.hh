// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/ESV2007.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/stationary-testcase.hh>
#include <dune/gdt/test/grids.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


// forward
template <class G, class R = double, int r = 1>
class ESV2007DdSubdomainsTestCase;


template <class E, class D, int d, class R, int r = 1>
class ESV2007Problem : public ProblemBase<E, D, d, R, r>
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class ESV2007Problem<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
    : public ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ScalarConstantFunctionType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> MatrixConstantFunctionType;
  typedef XT::Functions::ESV2007::Testcase1Force<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ForceType;

public:
  static const size_t default_integration_order = 2;

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration cfg;
    cfg["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    cfg["lower_left"] = "[-1 -1]";
    cfg["upper_right"] = "[1 1]";
    return cfg;
  }

  static XT::Common::Configuration default_boundary_info_cfg()
  {
    return XT::Grid::alldirichlet_boundaryinfo_default_config();
  }

  ESV2007Problem(const size_t integration_order = default_integration_order,
                 const XT::Common::Configuration& grd_cfg = default_grid_cfg(),
                 const XT::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(
          new ScalarConstantFunctionType(1, "diffusion_factor"),
          new MatrixConstantFunctionType(XT::Functions::internal::unit_matrix<RangeFieldImp, 2>(), "diffusion_tensor"),
          new ForceType(integration_order, "force"),
          new ScalarConstantFunctionType(0, "dirichlet"),
          new ScalarConstantFunctionType(0, "neumann"),
          grd_cfg,
          bnd_cfg)
  {
  }
}; // class ESV2007Problem< ..., 1 >


template <class G, class R = double, int r = 1>
class ESV2007TestCase
    : public Test::StationaryTestCase<G,
                                      LinearElliptic::ESV2007Problem<typename G::template Codim<0>::Entity,
                                                                     typename G::ctype,
                                                                     G::dimension,
                                                                     R,
                                                                     r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef XT::Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r> ExactSolutionType;

public:
  typedef LinearElliptic::ESV2007Problem<E, D, d, R, r> ProblemType;

private:
  typedef Test::StationaryTestCase<G, ProblemType> BaseType;

#if DXT_DISABLE_LARGE_TESTS

  template <class T, bool anything = true>
  struct Helper
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[4 4]";
      cfg["num_refinements"] = "0";
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
      cfg["num_elements"] = "[8 8]";
      cfg["num_refinements"] = "0";
      return cfg;
    }
  };

#if HAVE_DUNE_SPGRID
  template <class ct, int dim, template <int> class Ref, class Comm, bool anything>
  struct Helper<SPGrid<ct, dim, Ref, Comm>, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[8 8]";
      return cfg;
    }
  };
#endif

#if HAVE_DUNE_ALUGRID

  template <bool anything>
  struct Helper<AluConform2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[4 4]";
      cfg["num_refinements"] = "2";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluCube2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[8 8]";
      cfg["num_refinements"] = "1";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluSimplex2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[8 8]";
      cfg["num_refinements"] = "0";
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
  ESV2007TestCase(const size_t num_refs =
#if DXT_DISABLE_LARGE_TESTS
                      1
#else
                      3
#endif
                  )
    : BaseType(grid_cfg(), num_refs)
    , problem_()
    , exact_solution_()
  {
  }

  const ProblemType& problem() const override final
  {
    return problem_;
  }

  void print_header(std::ostream& out = DXTC_LOG_INFO_0) const override final
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase ESV2007: smooth data, homogeneous dirichlet          ||\n"
        << "||  (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)  ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
        << "||  diffusion = 1                                                 ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)              ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  exact solution = cos(1/2 pi x) cos(1/2 pi y)                  ||\n"
        << "|+================================================================+|\n"
        << "+==================================================================+" << std::endl;
  }

  bool provides_exact_solution() const override final
  {
    return true;
  }

  const ExactSolutionType& exact_solution() const override final
  {
    return exact_solution_;
  }

private:
  template <class G_, class R_, int r_>
  friend class ESV2007DdSubdomainsTestCase;

  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class ESV2007TestCase


template <class G, class R, int r>
class ESV2007DdSubdomainsTestCase
    : public Test::StationaryTestCase<G,
                                      LinearElliptic::ESV2007Problem<typename G::template Codim<0>::Entity,
                                                                     typename G::ctype,
                                                                     G::dimension,
                                                                     R,
                                                                     r>,
                                      XT::Grid::DD::SubdomainGrid<G>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef XT::Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r> ExactSolutionType;

public:
  typedef LinearElliptic::ESV2007Problem<E, D, d, R, r> ProblemType;

private:
  typedef Test::StationaryTestCase<G, ProblemType, XT::Grid::DD::SubdomainGrid<G>> BaseType;

  static XT::Common::Configuration grid_cfg()
  {
    auto cfg = ESV2007TestCase<G, R, r>::grid_cfg();
    cfg["type"] = XT::Grid::cube_dd_subdomains_gridprovider_id();
    cfg["num_partitions"] = "[1 1 1 1]";
    return cfg;
  }

public:
  ESV2007DdSubdomainsTestCase(const size_t num_refs =
#if DXT_DISABLE_LARGE_TESTS
                                  1
#else
                                  3
#endif
                              )
    : BaseType(grid_cfg(), num_refs)
    , problem_()
    , exact_solution_()
  {
  }

  const ProblemType& problem() const override final
  {
    return problem_;
  }

  void print_header(std::ostream& out = DXTC_LOG_INFO_0) const override final
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase ESV2007: smooth data, homogeneous dirichlet          ||\n"
        << "||    DD::SubdomainGrid variant                                   ||\n"
        << "||  (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)  ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
        << "||  diffusion = 1                                                 ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)              ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  exact solution = cos(1/2 pi x) cos(1/2 pi y)                  ||\n"
        << "|+================================================================+|\n"
        << "+==================================================================+" << std::endl;
  }

  bool provides_exact_solution() const override final
  {
    return true;
  }

  const ExactSolutionType& exact_solution() const override final
  {
    return exact_solution_;
  }

private:
  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class ESV2007DdSubdomainsTestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
