// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ER2007_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ER2007_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/expression.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/stationary-testcase.hh>
#include <dune/gdt/test/grids.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template <class E, class D, int d, class R, int r = 1>
class ER2007Problem : public ProblemBase<E, D, d, R, r>
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class ER2007Problem<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
    : public ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ScalarConstantFunctionType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> MatrixConstantFunctionType;
  typedef XT::Functions::ExpressionFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ExpressionFunctionType;

public:
  static const size_t default_integration_order = 3;

  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration cfg;
    cfg["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    cfg["lower_left"] = "[0 0]";
    cfg["upper_right"] = "[1 1]";
    return cfg;
  }

  static XT::Common::Configuration default_boundary_info_cfg()
  {
    return XT::Grid::alldirichlet_boundaryinfo_default_config();
  }

  ER2007Problem(const size_t integration_order = default_integration_order,
                const XT::Common::Configuration& grd_cfg = default_grid_cfg(),
                const XT::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(
          new ScalarConstantFunctionType(1, "diffusion_factor"),
          new MatrixConstantFunctionType(XT::Functions::internal::unit_matrix<RangeFieldImp, 2>(), "diffusion_tensor"),
          new ExpressionFunctionType("x", "64.0*pi*pi*(cos(8.0*pi*x[0])+cos(8.0*pi*x[1]))", integration_order, "force"),
          new ExpressionFunctionType("x", "cos(8.0*pi*x[0])+cos(8.0*pi*x[1])", integration_order, "dirichlet"),
          new ScalarConstantFunctionType(0, "neumann"),
          grd_cfg,
          bnd_cfg)
  {
  }
}; // class ER2007Problem< ..., 1 >


template <class G, class R = double, int r = 1>
class ER2007TestCase
    : public Test::StationaryTestCase<G,
                                      LinearElliptic::ER2007Problem<typename G::template Codim<0>::Entity,
                                                                    typename G::ctype,
                                                                    G::dimension,
                                                                    R,
                                                                    r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef XT::Functions::ExpressionFunction<E, D, d, R, r> ExactSolutionType;

public:
  typedef LinearElliptic::ER2007Problem<E, D, d, R, r> ProblemType;

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
      cfg["num_elements"] = "[8 8]";
      return cfg;
    }
  };

#if HAVE_ALUGRID
  template <bool anything>
  struct Helper<AluConform2dGridType, anything>
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
      return cfg;
    }
  };

#endif // HAVE_ALUGRID
#endif // DXT_DISABLE_LARGE_TESTS

  static XT::Common::Configuration grid_cfg()
  {
    auto cfg = ProblemType::default_grid_cfg();
    cfg = Helper<typename std::decay<G>::type>::value(cfg);
    return cfg;
  }

public:
  using typename BaseType::GridType;

  ER2007TestCase(const size_t num_refs =
#if DXT_DISABLE_LARGE_TESTS
                     1
#else
                     2
#endif
                 )
    : BaseType(XT::Grid::make_cube_grid<GridType>(grid_cfg()).grid_ptr(), num_refs)
    , problem_()
    , exact_solution_("x",
                      "cos(8.0*pi*x[0])+cos(8.0*pi*x[1])",
                      ProblemType::default_integration_order,
                      "exact solution",
                      {"-8.0*pi*sin(8.0*pi*x[0])", "-8.0*pi*sin(8.0*pi*x[1])"})
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+============================================================+\n"
        << "|+==========================================================+|\n"
        << "||  Testcase ER2007: smooth data, nonhomogeneous dirichlet  ||\n"
        << "||  (see page 858 in Epshteyn, Riviere, 2007)               ||\n"
        << "|+----------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                                ||\n"
        << "||  diffusion = 1                                           ||\n"
        << "||  force     = 64 pi^2 (cos(8 pi x) + cos(8 pi y))         ||\n"
        << "||  dirichlet = cos(8 pi x) + cos(8 pi y)                   ||\n"
        << "||  exact solution = cos(8 pi x) + cos(8 pi y)              ||\n"
        << "|+==========================================================+|\n"
        << "+============================================================+" << std::endl;
  }

  virtual bool provides_exact_solution() const override final
  {
    return true;
  }

  virtual const ExactSolutionType& exact_solution() const override final
  {
    return exact_solution_;
  }

private:
  const ProblemType problem_;
  const ExactSolutionType exact_solution_;
}; // class ER2007TestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ER2007_HH
