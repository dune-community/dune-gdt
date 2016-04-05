// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/ESV2007.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/test/stationary-testcase.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template< class E, class D, int d, class R, int r = 1 >
class ESV2007Problem
  : public ProblemBase< E, D, d, R, r >
{
  static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class ESV2007Problem< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >                   BaseType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > MatrixConstantFunctionType;
  typedef Stuff::Functions::ESV2007::Testcase1Force< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 > ForceType;

public:
  static const size_t default_integration_order = 2;

  static Stuff::Common::Configuration default_grid_cfg()
  {
    Stuff::Common::Configuration cfg;
    cfg["type"] = Stuff::Grid::Providers::Configs::Cube_default()["type"];
    cfg["lower_left"]  = "[-1 -1]";
    cfg["upper_right"] = "[1 1]";
    return cfg;
  }

  static Stuff::Common::Configuration default_boundary_info_cfg()
  {
    return Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();
  }

  ESV2007Problem(const size_t integration_order = default_integration_order,
                 const Stuff::Common::Configuration& grd_cfg = default_grid_cfg(),
                 const Stuff::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(new ScalarConstantFunctionType(1, "diffusion_factor"),
               new MatrixConstantFunctionType(Stuff::Functions::internal::unit_matrix< RangeFieldImp, 2 >(), "diffusion_tensor"),
               new ForceType(integration_order,  "force"),
               new ScalarConstantFunctionType(0, "dirichlet"),
               new ScalarConstantFunctionType(0, "neumann"),
               grd_cfg,
               bnd_cfg)
  {}
}; // class ESV2007Problem< ..., 1 >


template< class G, class R = double, int r = 1 >
class ESV2007TestCase
  : public Test::StationaryTestCase< G,
                                     LinearElliptic::ESV2007Problem< typename G::template Codim< 0 >::Entity,
                                                                     typename G::ctype, G::dimension, R, r > >
{
  typedef typename G::template Codim< 0 >::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef Stuff::Functions::ESV2007::Testcase1ExactSolution< E, D, d, R, r > ExactSolutionType;
public:
  typedef LinearElliptic::ESV2007Problem< E, D, d, R, r > ProblemType;
private:
  typedef Test::StationaryTestCase< G, ProblemType > BaseType;

  template< class T, bool anything = true >
  struct Helper
  {
    static_assert(AlwaysFalse< T >::value, "Please add a configuration for this grid type!");
    static Stuff::Common::Configuration value(Stuff::Common::Configuration cfg) { return cfg; }
  };

  template< bool anything >
  struct Helper< SGrid< 2, 2 >, anything >
  {
    static Stuff::Common::Configuration value(Stuff::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[8 8]";
      return cfg;
    }
  };

#if HAVE_ALUGRID
  template< bool anything >
  struct Helper< ALUGrid< 2, 2, simplex, conforming >, anything >
  {
    static Stuff::Common::Configuration value(Stuff::Common::Configuration cfg)
    {
      cfg["num_elements"]    = "[4 4]";
      cfg["num_refinements"] = "2";
      return cfg;
    }
  };

  template< bool anything >
  struct Helper< ALUGrid< 2, 2, simplex, nonconforming >, anything >
  {
    static Stuff::Common::Configuration value(Stuff::Common::Configuration cfg)
    {
      cfg["num_elements"] = "[8 8]";
      return cfg;
    }
  };
#endif // HAVE_ALUGRID

  static Stuff::Common::Configuration grid_cfg()
  {
    auto cfg = ProblemType::default_grid_cfg();
    cfg = Helper< typename std::decay< G >::type >::value(cfg);
    return cfg;
  }
public:
  using typename BaseType::GridType;

  ESV2007TestCase(const size_t num_refs = 3)
    : BaseType(Stuff::Grid::Providers::Cube< G >::create(grid_cfg())->grid_ptr(), num_refs)
    , problem_()
    , exact_solution_()
  {}

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
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
}; // class ESV2007TestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ESV2007_HH
