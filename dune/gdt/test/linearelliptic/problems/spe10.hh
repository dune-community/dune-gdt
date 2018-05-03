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

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/stationary-testcase.hh>
#include <dune/gdt/test/grids.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template <class E, class D, int d, class R, int r = 1>
class Spe10Model1Problem : public ProblemBase<E, D, d, R, r>
{
  static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Spe10Model1Problem<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
    : public ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef ProblemBase<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ScalarConstantFunctionType;
  typedef XT::Functions::ConstantFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> MatrixConstantFunctionType;
  typedef XT::Functions::IndicatorFunction<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> IndicatorFunctionType;
  typedef XT::Functions::Spe10::Model1Function<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> Spe10FunctionType;

public:
  static XT::Common::Configuration default_grid_cfg()
  {
    XT::Common::Configuration cfg;
    cfg["type"] = XT::Grid::cube_gridprovider_default_config()["type"];
    cfg["lower_left"] = "[0 0]";
    cfg["upper_right"] = "[5 1]";
    return cfg;
  }

  static XT::Common::Configuration default_boundary_info_cfg()
  {
    return XT::Grid::alldirichlet_boundaryinfo_default_config();
  }

  Spe10Model1Problem(const XT::Common::Configuration& grd_cfg = default_grid_cfg(),
                     const XT::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(
          new Spe10FunctionType(XT::Functions::Spe10::internal::model1_filename,
                                grd_cfg.get<typename Spe10FunctionType::DomainType>("lower_left"),
                                grd_cfg.get<typename Spe10FunctionType::DomainType>("upper_right"),
                                XT::Functions::Spe10::internal::model1_min_value,
                                XT::Functions::Spe10::internal::model1_max_value,
                                "diffusion_factor"),
          new MatrixConstantFunctionType(XT::Functions::internal::unit_matrix<RangeFieldImp, 2>(), "diffusion_tensor"),
          new IndicatorFunctionType({{{{0.95, 0.30}, {1.10, 0.45}}, 2000},
                                     {{{3.00, 0.75}, {3.15, 0.90}}, -1000},
                                     {{{4.25, 0.25}, {4.40, 0.40}}, -1000}},
                                    "force"),
          new ScalarConstantFunctionType(0, "dirichlet"),
          new ScalarConstantFunctionType(0, "neumann"),
          grd_cfg,
          bnd_cfg)
  {
  }
}; // class Spe10Model1Problem< ..., 1 >


template <class G, class R = double, int r = 1>
class Spe10Model1TestCase
    : public Test::StationaryTestCase<G,
                                      LinearElliptic::Spe10Model1Problem<typename G::template Codim<0>::Entity,
                                                                         typename G::ctype,
                                                                         G::dimension,
                                                                         R,
                                                                         r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  typedef LinearElliptic::Spe10Model1Problem<E, D, d, R, r> ProblemType;

private:
  typedef Test::StationaryTestCase<G, ProblemType> BaseType;

#if DXT_DISABLE_LARGE_TESTS

  template <class T, bool anything = true>
  struct Helper
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[20 4]";
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
      cfg["num_elements"] = "[100 20]";
      return cfg;
    }
  };

#if HAVE_DUNE_ALUGRID

  template <bool anything>
  struct Helper<AluConform2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[100 20]";
      cfg["num_refinements"] = "1";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluCube2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[100 20]";
      cfg["num_refinements"] = "1";
      return cfg;
    }
  };

  template <bool anything>
  struct Helper<AluSimplex2dGridType, anything>
  {
    static XT::Common::Configuration value(XT::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[100 20]";
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

  Spe10Model1TestCase(const size_t num_refs = 1)
    : BaseType(grid_cfg(), num_refs)
    , problem_()
  {
  }

  Spe10Model1TestCase(XT::Common::Configuration cfg)
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
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: SPE10, Model1                               ||\n"
        << "||  (see http://www.spe.org/web/csp/datasets/set01.htm)   ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 5] x [0, 1]                              ||\n"
        << "||  diffusion: spe10 model 1 scalar data                  ||\n"
        << "||         |  2000 in [0.95, 1.10] x [0.30, 0.45]         ||\n"
        << "||  force: | -1000 in [3.00, 3.15] x [0.75, 0.90]         ||\n"
        << "||         | -1000 in [4.25, 4.40] x [0.25, 0.40]         ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class Spe10Model1TestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH
