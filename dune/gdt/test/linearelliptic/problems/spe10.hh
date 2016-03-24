// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH

#include <dune/grid/alugrid.hh>

#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/indicator.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/test/stationary-eocstudy.hh>

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
  typedef Stuff::Functions::Constant<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> ScalarConstantFunctionType;
  typedef Stuff::Functions::Constant<EntityImp, DomainFieldImp, 2, RangeFieldImp, 2, 2> MatrixConstantFunctionType;
  typedef Stuff::Functions::Indicator<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> IndicatorFunctionType;
  typedef Stuff::Functions::Spe10::Model1<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> Spe10FunctionType;

public:
  static Stuff::Common::Configuration default_grid_cfg()
  {
    auto cfg               = Stuff::Grid::Providers::Configs::Cube_default();
    cfg["lower_left"]      = "[0 0]";
    cfg["upper_right"]     = "[5 1]";
    cfg["num_elements"]    = "[100 20]";
    cfg["num_refinements"] = "1";
    return cfg;
  }

  static Stuff::Common::Configuration default_boundary_info_cfg()
  {
    return Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config();
  }

  Spe10Model1Problem(const Stuff::Common::Configuration& grd_cfg = default_grid_cfg(),
                     const Stuff::Common::Configuration& bnd_cfg = default_boundary_info_cfg())
    : BaseType(new Spe10FunctionType(Stuff::Functions::Spe10::internal::model1_filename,
                                     grd_cfg.get<typename Spe10FunctionType::DomainType>("lower_left"),
                                     grd_cfg.get<typename Spe10FunctionType::DomainType>("upper_right"),
                                     Stuff::Functions::Spe10::internal::model1_min_value,
                                     Stuff::Functions::Spe10::internal::model1_max_value, "diffusion_factor"),
               new MatrixConstantFunctionType(Stuff::Functions::internal::unit_matrix<RangeFieldImp, 2>(),
                                              "diffusion_tensor"),
               new IndicatorFunctionType({{{{0.95, 0.30}, {1.10, 0.45}}, 2000},
                                          {{{3.00, 0.75}, {3.15, 0.90}}, -1000},
                                          {{{4.25, 0.25}, {4.40, 0.40}}, -1000}},
                                         "force"),
               new ScalarConstantFunctionType(0, "dirichlet"), new ScalarConstantFunctionType(0, "neumann"), grd_cfg,
               bnd_cfg)
  {
  }
}; // class Spe10Model1Problem< ..., 1 >


template <class G, class R = double, int r = 1>
class Spe10Model1TestCase
    : public Test::StationaryTestCase<G, LinearElliptic::Spe10Model1Problem<typename G::template Codim<0>::Entity,
                                                                            typename G::ctype, G::dimension, R, r>>
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;

public:
  typedef LinearElliptic::Spe10Model1Problem<E, D, d, R, r> ProblemType;

private:
  typedef Test::StationaryTestCase<G, ProblemType> BaseType;

public:
  using typename BaseType::GridType;

  Spe10Model1TestCase(const size_t num_refs = 1)
    : BaseType(create_initial_grid(), num_refs)
    , problem_()
  {
  }

  virtual const ProblemType& problem() const override final
  {
    return problem_;
  }

  virtual void print_header(std::ostream& out = std::cout) const override final
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: SPE10, Model1                               ||\n"
        << "||  (see http://www.spe.org/web/csp/datasets/set01.htm)   ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 5] x [0, 1]                              ||\n"
        << "||  diffusion: spe10 model 1 scalar data                  ||\n"
        << "||         |  2000 in [0.55, 0.70] x [0.70, 0.85]         ||\n"
        << "||  force: | -1000 in [3.00, 3.15] x [0.77, 0.90]         ||\n"
        << "||         | -1000 in [4.30, 4.45] x [0.50, 0.65]         ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  }

private:
  static std::shared_ptr<GridType> create_initial_grid()
  {
    auto grid_cfg = ProblemType::default_grid_cfg();
#if HAVE_ALUGRID
    if (std::is_same<GridType, ALUGrid<2, 2, simplex, conforming>>::value)
      grid_cfg["num_refinements"] = "1";
#endif // HAVE_ALUGRID
    return Stuff::Grid::Providers::Cube<G>::create(grid_cfg)->grid_ptr();
  } // ... create_initial_grid(...)

  const ProblemType problem_;
}; // class Spe10Model1TestCase


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_PROBLEMS_SPE10_HH
