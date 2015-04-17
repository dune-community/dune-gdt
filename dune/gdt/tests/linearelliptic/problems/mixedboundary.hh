// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_MIXEDBOUNDARY_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_MIXEDBOUNDARY_HH

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
class MixedBoundaryProblem
  : public ProblemBase< E, D, d, R, r >
{
  static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!");
};


template< class EntityImp, class DomainFieldImp, class RangeFieldImp >
class MixedBoundaryProblem< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
  : public ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >
{
  typedef ProblemBase< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >            BaseType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >    ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant< EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2 > MatrixConstantFunctionType;
  typedef Stuff::Functions::Expression< EntityImp, DomainFieldImp, 2, RangeFieldImp, 1 >  ExpressionFunctionType;

public:
  static const size_t default_integration_order = 2;

  static Stuff::Common::Configuration default_grid_cfg()
  {
    auto cfg = Stuff::Grid::Providers::Configs::Cube_default();
    cfg["num_elements"]    = "[2 2]";
    cfg["num_refinements"] = "1";
    return cfg;
  }

  static Stuff::Common::Configuration default_boundary_info_cfg()
  {
    return Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();
  }

  MixedBoundaryProblem(const size_t integration_order = default_integration_order,
                 const Stuff::Common::Configuration& grd_cfg = default_grid_cfg(),
                 const Stuff::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(new ScalarConstantFunctionType(1, "diffusion_factor"),
               new MatrixConstantFunctionType(Stuff::Functions::internal::unit_matrix< RangeFieldImp, 2 >(), "diffusion_tensor"),
               new ScalarConstantFunctionType(1,  "force"),
               new ExpressionFunctionType("x", "0.25 * x[0] * x[1]", integration_order, "dirichlet"),
               new ScalarConstantFunctionType(0.1, "neumann"),
               grd_cfg,
               bnd_cfg)
  {}
}; // class MixedBoundaryProblem< ..., 1 >


template< class G, class R = double, int r = 1 >
class MixedBoundaryTestCase
  : public Tests::StationaryTestCase< G,
                                      LinearElliptic::MixedBoundaryProblem< typename G::template Codim< 0 >::Entity,
                                                                      typename G::ctype, G::dimension, R, r > >
{
  typedef typename G::template Codim< 0 >::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
public:
  typedef LinearElliptic::MixedBoundaryProblem< E, D, d, R, r > ProblemType;
private:
  typedef Tests::StationaryTestCase< G, ProblemType > BaseType;
public:
  using typename BaseType::GridType;

  MixedBoundaryTestCase(const size_t num_refs = 3)
    : BaseType(Stuff::Grid::Providers::Cube< G >::create(ProblemType::default_grid_cfg())->grid_ptr(), num_refs)
    , problem_()
  {}

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: mixed boundary types                         ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                              ||\n"
        << "||  diffusion = 1                                         ||\n"
        << "||  force     = 1                                         ||\n"
        << "||  neumann   = 0.1       on the right side               ||\n"
        << "||  dirichlet = 1/4 x y   everywhere else                 ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  }

private:
  const ProblemType problem_;
}; // class MixedBoundaryTestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_MIXEDBOUNDARY_HH
