// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ER2007_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_ER2007_HH

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/tests/stationary-eocstudy.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


template< class E, class D, int d, class R, int r = 1 >
class ER2007Problem
  : public ProblemBase< E, D, d, R, r >
{
  static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class ER2007Problem< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >                   BaseType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > MatrixConstantFunctionType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >  ExpressionFunctionType;

public:
  static const size_t default_integration_order = 3;

  static Stuff::Common::Configuration default_grid_cfg()
  {
    auto cfg = Stuff::Grid::Providers::Configs::Cube_default();
    cfg["num_elements"]    = "[16 16]";
    cfg["num_refinements"] = "1";
    return cfg;
  }

  static Stuff::Common::Configuration default_boundary_info_cfg()
  {
    return Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();
  }

  ER2007Problem(const size_t integration_order = default_integration_order,
                 const Stuff::Common::Configuration& grd_cfg = default_grid_cfg(),
                 const Stuff::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(new ScalarConstantFunctionType(1, "diffusion_factor"),
               new MatrixConstantFunctionType(Stuff::Functions::internal::unit_matrix< RangeFieldImp, 2 >(),
                                              "diffusion_tensor"),
               new ExpressionFunctionType("x",
                                          "64.0*pi*pi*(cos(8.0*pi*x[0])+cos(8.0*pi*x[1]))",
                                          integration_order,
                                          "force"),
               new ExpressionFunctionType("x",
                                          "cos(8.0*pi*x[0])+cos(8.0*pi*x[1])",
                                          integration_order,
                                          "dirichlet"),
               new ScalarConstantFunctionType(0, "neumann"),
               grd_cfg,
               bnd_cfg)
  {}
}; // class ER2007Problem< ..., 1 >


template< class G, class R = double, int r = 1 >
class ER2007TestCase
  : public Tests::StationaryTestCase< G,
                                      LinearElliptic::ER2007Problem< typename G::template Codim< 0 >::Entity,
                                                                      typename G::ctype, G::dimension, R, r > >
{
  typedef typename G::template Codim< 0 >::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef Stuff::Functions::Expression< E, D, d, R, r > ExactSolutionType;
public:
  typedef LinearElliptic::ER2007Problem< E, D, d, R, r > ProblemType;
private:
  typedef Tests::StationaryTestCase< G, ProblemType > BaseType;
public:
  using typename BaseType::GridType;

  ER2007TestCase(const size_t num_refs = 1)
    : BaseType(Stuff::Grid::Providers::Cube< G >::create(ProblemType::default_grid_cfg())->grid_ptr(), num_refs)
    , problem_()
    , exact_solution_("x",
                      "cos(8.0*pi*x[0])+cos(8.0*pi*x[1])",
                      ProblemType::default_integration_order,
                      "exact solution",
                      {{"-8.0*pi*sin(8.0*pi*x[0])",
                        "-8.0*pi*sin(8.0*pi*x[1])"}})
  {}

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
