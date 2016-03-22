// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_HH
#define DUNE_GDT_TEST_ELLIPTIC_HH

#define DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS 1

#ifdef HAVE_FASP
# undef HAVE_FASP
#endif

#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/dgfparser.hh>
#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/gdt/spaces/tools.hh>

std::ostream&
  DUNE_DEPRECATED_MSG("Use DSC_LOG_INFO instead and toggle DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING")
              test_out = std::cout;

namespace EllipticTestCase {


template< class GridType >
class Base
{
  typedef DSG::ProviderInterface< GridType > GridProviderType;
public:
#if HAVE_DUNE_FEM
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< GridType, false >::Type  GridPartType;
#endif
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< GridType >::Type         GridViewType;
  typedef typename GridType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridType::ctype                        DomainFieldType;
  static const size_t                                     dimDomain = GridType::dimension;

  Base(std::unique_ptr< GridProviderType >&& grd_prv, size_t num_refinements)
    : grid_provider_(std::move(grd_prv))
  {
    levels_.push_back(grid().maxLevel());
#if HAVE_DUNE_FEM
    level_grid_parts_.emplace_back(grid(), grid().maxLevel());
#endif
    level_grid_views_.emplace_back(grid().levelGridView(grid().maxLevel()));
    for (size_t rr = 0; rr < num_refinements; ++rr) {
      grid().globalRefine(Dune::DGFGridInfo< GridType >::refineStepsForHalf());
      levels_.push_back(grid().maxLevel());
#if HAVE_DUNE_FEM
      level_grid_parts_.emplace_back(grid(), grid().maxLevel());
#endif
      level_grid_views_.emplace_back(grid().levelGridView(grid().maxLevel()));
    }
    grid().globalRefine(Dune::DGFGridInfo< GridType >::refineStepsForHalf());
#if HAVE_DUNE_FEM
    reference_grid_part_ = DSC::make_unique< GridPartType >(grid(), grid().maxLevel());
#endif
    reference_grid_view_ = DSC::make_unique< GridViewType >(grid().levelGridView(grid().maxLevel()));
  }

  size_t num_levels() const
  {
    return levels_.size();
  }

#if HAVE_DUNE_FEM
  GridPartType level_grid_part(const size_t level) const
  {
    assert(level < levels_.size());
    assert(level < level_grid_parts_.size());
    assert(ssize_t(levels_[level]) < grid().maxLevel());
    return level_grid_parts_[level];
  }
#endif // HAVE_DUNE_FEM

  GridViewType level_grid_view(const size_t level) const
  {
    assert(level < levels_.size());
    assert(level < level_grid_views_.size());
    assert(ssize_t(levels_[level]) < grid().maxLevel());
    return level_grid_views_[level];
  }

#if HAVE_DUNE_FEM
  GridPartType reference_grid_part() const
  {
    return *reference_grid_part_;
  }
#endif // HAVE_DUNE_FEM

  GridViewType reference_grid_view() const
  {
    return *reference_grid_view_;
  }

private:
  GridType& grid() { return grid_provider_->grid(); }
  const GridType& grid() const { return grid_provider_->grid(); }


  std::unique_ptr< GridProviderType > grid_provider_;
  std::vector< size_t > levels_;
#if HAVE_DUNE_FEM
  std::vector< GridPartType > level_grid_parts_;
  std::unique_ptr< GridPartType > reference_grid_part_;
#endif
  std::vector< GridViewType > level_grid_views_;
  std::unique_ptr< GridViewType > reference_grid_view_;
}; // class Base


template< class GridType >
class ER07
  : public Base< GridType >
{
  typedef DSG::Providers::Cube< GridType > GridProviderType;
  typedef Base< GridType > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef double            RangeFieldType;
  static const size_t       dimRange = 1;
  typedef DSG::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;

  typedef Dune::Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ExpressionFunctionType;
  typedef ConstantFunctionType    DiffusionType;
  typedef ExpressionFunctionType  ForceType;
  typedef ExpressionFunctionType  DirichletType;
  typedef ConstantFunctionType    NeumannType;
  typedef ExpressionFunctionType  ExactSolutionType;

  ER07(const size_t num_refinements = 2)
    : BaseType(create_initial_grid(), num_refinements)
    , boundary_info_()
    , diffusion_(1)
    , force_("x", "64.0 * pi * pi * (cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1]))", 3)
    , dirichlet_("x", "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])", 3)
    , neumann_(0)
    , exact_solution_("x",
                      "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])",
                      3,
                      "exact solution",
                      {"-8.0 * pi * sin(8.0 * pi * x[0])",
                       "-8.0 * pi * sin(8.0 * pi * x[1])"})
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase ER07: smooth data, nonhomogeneous dirichlet  ||\n"
        << "||  (see page 858 in Epshteyn, Riviere, 2007)             ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                              ||\n"
        << "||  diffusion = 1                                         ||\n"
        << "||  force     = 64 pi^2 (cos(8 pi x) + cos(8 pi y))       ||\n"
        << "||  dirichlet = cos(8 pi x) + cos(8 pi y)                 ||\n"
        << "||  exact solution = cos(8 pi x) + cos(8 pi y)            ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  } // ... print_header(...)

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const DiffusionType& diffusion() const
  {
    return diffusion_;
  }

  const ForceType& force() const
  {
    return force_;
  }

  const DirichletType& dirichlet() const
  {
    return dirichlet_;
  }

  const NeumannType& neumann() const
  {
    return neumann_;
  }

  bool provides_exact_solution() const
  {
    return true;
  }

  const ExactSolutionType& exact_solution() const
  {
    return exact_solution_;
  }

private:
  static std::unique_ptr< GridProviderType > create_initial_grid()
  {
    auto grid_provider = DSC::make_unique< GridProviderType >(0, 1, 16);
    auto& grid = grid_provider->grid();
    grid.globalRefine(1);
    return grid_provider;
  } // ... create_initial_grid(...)

  const BoundaryInfoType boundary_info_;
  const DiffusionType diffusion_;
  const ForceType force_;
  const DirichletType dirichlet_;
  const NeumannType neumann_;
  const ExactSolutionType exact_solution_;
}; // class ER07


template< class GridType >
class ESV07
  : public Base< GridType >
{
  typedef DSG::Providers::Cube< GridType > GridProviderType;
  typedef Base< GridType > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef double            RangeFieldType;
  static const size_t       dimRange = 1;
  typedef DSG::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;

  typedef Dune::Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ExpressionFunctionType;
  typedef ConstantFunctionType    DiffusionType;
  typedef ExpressionFunctionType  ForceType;
  typedef ConstantFunctionType    DirichletType;
  typedef ConstantFunctionType    NeumannType;
  typedef ExpressionFunctionType  ExactSolutionType;

  ESV07(const size_t num_refinements = 3)
    : BaseType(create_initial_grid(), num_refinements)
    , boundary_info_()
    , diffusion_(1)
    , force_("x", "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])", 2)
    , dirichlet_(0)
    , neumann_(0)
    , exact_solution_("x",
                      "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                      2,
                      "exact solution",
                      {"-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                       "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"})
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==================================================================+\n"
        << "|+================================================================+|\n"
        << "||  Testcase ESV07: smooth data, homogeneous dirichlet            ||\n"
        << "||  (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)  ||\n"
        << "|+----------------------------------------------------------------+|\n"
        << "||  domain = [-1, 1] x [-1, 1]                                    ||\n"
        << "||  diffusion = 1                                                 ||\n"
        << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)              ||\n"
        << "||  dirichlet = 0                                                 ||\n"
        << "||  exact solution = cos(1/2 pi x) cos(1/2 pi y)                  ||\n"
        << "|+================================================================+|\n"
        << "+==================================================================+" << std::endl;
  } // ... print_header(...)

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const DiffusionType& diffusion() const
  {
    return diffusion_;
  }

  const ForceType& force() const
  {
    return force_;
  }

  const DirichletType& dirichlet() const
  {
    return dirichlet_;
  }

  const NeumannType& neumann() const
  {
    return neumann_;
  }

  bool provides_exact_solution() const
  {
    return true;
  }

  const ExactSolutionType& exact_solution() const
  {
    return exact_solution_;
  }

private:
  static std::unique_ptr< GridProviderType > create_initial_grid()
  {
    if (std::is_same< GridType, Dune::SGrid< 2, 2 > >::value) {
      return DSC::make_unique<GridProviderType>(-1, 1, 8);
#if HAVE_ALUGRID
    } else if (std::is_same< GridType, Dune::ALUConformGrid< 2, 2 > >::value
        || std::is_same< GridType, Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >::value) {
      auto grid_provider = DSC::make_unique<DSG::Providers::Cube< GridType >>(-1, 1, 4);
      grid_provider->grid().globalRefine(2);
      return grid_provider;
#endif // HAVE_ALUGRID
    } else {
      auto grid_provider = DSC::make_unique<GridProviderType>(-1, 1, 4);
      grid_provider->grid().globalRefine(1);
      return grid_provider;
    }
  } // ... create_initial_grid(...)

  const BoundaryInfoType boundary_info_;
  const DiffusionType diffusion_;
  const ForceType force_;
  const DirichletType dirichlet_;
  const NeumannType neumann_;
  const ExactSolutionType exact_solution_;
}; // class ESV07


template< class GridType >
class LocalThermalBlock
  : public Base< GridType >
{
  typedef DSG::Providers::Cube< GridType > GridProviderType;
  typedef Base< GridType > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef double            RangeFieldType;
  static const size_t       dimRange = 1;
  typedef DSG::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;

  typedef Dune::Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Dune::Stuff::Functions::Checkerboard
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    CheckerboardFunctionType;
  typedef CheckerboardFunctionType  DiffusionType;
  typedef ConstantFunctionType      ForceType;
  typedef ConstantFunctionType      DirichletType;
  typedef ConstantFunctionType      NeumannType;
  typedef ConstantFunctionType      ExactSolutionType;

  LocalThermalBlock(const size_t num_refinements = 3)
    : BaseType(create_initial_grid(), num_refinements)
    , boundary_info_()
    , diffusion_({0.0, 0.0},
                 {1.0, 1.0},
                 {6, 6},
                 {1.0,  1.0, 1.0, 0.1, 0.1, 0.1,
                  1.0, 0.01, 1.0, 0.1, 0.1, 0.1,
                  1.0,  1.0, 1.0, 0.1, 0.1, 0.1,
                  1.0,  1.0, 1.0, 0.1, 0.1, 0.1,
                  1.0, 0.01, 1.0, 0.1, 0.1, 0.1,
                  1.0,  1.0, 1.0, 0.1, 0.1, 0.1})
    , force_(1)
    , dirichlet_(0)
    , neumann_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+======================================================================+\n"
        << "|+====================================================================+|\n"
        << "||  Testcase: local thermal block problem                             ||\n"
        << "||  (see http://wwwmath.uni-muenster.de/num/publications/2013/AO13/)  ||\n"
        << "|+--------------------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                                          ||\n"
        << "||  diffusion:  see page 3 (mu_test)                                  ||\n"
        << "||  force     = 1                                                     ||\n"
        << "||  dirichlet = 0                                                     ||\n"
        << "||  reference solution: discrete solution on finest grid              ||\n"
        << "|+====================================================================+|\n"
        << "+======================================================================+" << std::endl;
  } // ... print_header(...)

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const DiffusionType& diffusion() const
  {
    return diffusion_;
  }

  const ForceType& force() const
  {
    return force_;
  }

  const DirichletType& dirichlet() const
  {
    return dirichlet_;
  }

  const NeumannType& neumann() const
  {
    return neumann_;
  }

  bool provides_exact_solution() const
  {
    return false;
  }

  const ExactSolutionType& exact_solution() const
  {
    DUNE_THROW(Dune::InvalidStateException, "provides_exact_solution() == false");
    return neumann_;
  }

private:
  static std::unique_ptr< GridProviderType > create_initial_grid()
  {
    auto grid_provider = DSC::make_unique< GridProviderType >(0, 1, 6);
    grid_provider->grid().globalRefine(1);
    return grid_provider;
  } // ... create_initial_grid(...)

  const BoundaryInfoType boundary_info_;
  const DiffusionType diffusion_;
  const ForceType force_;
  const DirichletType dirichlet_;
  const NeumannType neumann_;
}; // class LocalThermalBlock


template< class GridType >
class MixedBoundaryTypes
  : public Base< GridType >
{
  typedef DSG::Providers::Cube< GridType > GridProviderType;
  typedef Base< GridType > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef double            RangeFieldType;
  static const size_t       dimRange = 1;
  typedef DSG::BoundaryInfos::NormalBased< typename GridViewType::Intersection > BoundaryInfoType;

  typedef Dune::Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ExpressionFunctionType;
  typedef ConstantFunctionType    DiffusionType;
  typedef ConstantFunctionType    ForceType;
  typedef ExpressionFunctionType  DirichletType;
  typedef ConstantFunctionType    NeumannType;
  typedef ConstantFunctionType    ExactSolutionType;

  MixedBoundaryTypes(const size_t num_refinements = 3)
    : BaseType(create_initial_grid(), num_refinements)
    , boundary_info_(create_boundary_info())
    , diffusion_(1)
    , force_(1)
    , dirichlet_("x", "0.25 * x[0] * x[1]", 2)
    , neumann_(0.1)
    , exact_solution_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase mixed boundary types                         ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 1] x [0, 1]                              ||\n"
        << "||  diffusion = 1                                         ||\n"
        << "||  force     = 1                                         ||\n"
        << "||  neumann   = 0.1       on the right side               ||\n"
        << "||  dirichlet = 1/4 x y   everywhere else                 ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  } // ... print_header(...)

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const DiffusionType& diffusion() const
  {
    return diffusion_;
  }

  const ForceType& force() const
  {
    return force_;
  }

  const DirichletType& dirichlet() const
  {
    return dirichlet_;
  }

  const NeumannType& neumann() const
  {
    return neumann_;
  }

  bool provides_exact_solution() const
  {
    return false;
  }

  const ExactSolutionType& exact_solution() const
  {
    DUNE_THROW(Dune::InvalidStateException, "provides_exact_solution() == false");
    return neumann_;
  }

private:
  static std::unique_ptr< GridProviderType > create_initial_grid()
  {
    auto grid_provider = DSC::make_unique< GridProviderType >(0, 1, 2);
    grid_provider->grid().globalRefine(1);
    return grid_provider;
  } // ... create_initial_grid(...)

  static BoundaryInfoType create_boundary_info()
  {
    typename ConstantFunctionType::DomainType neumann_normal(0.0);
    neumann_normal[0] = 1.0;
    return BoundaryInfoType(true, {}, {neumann_normal});
  }

  const BoundaryInfoType boundary_info_;
  const DiffusionType diffusion_;
  const ForceType force_;
  const DirichletType dirichlet_;
  const NeumannType neumann_;
  const ExactSolutionType exact_solution_;
}; // class MixedBoundaryTypes


template< class GridType >
class Spe10Model1
  : public Base< GridType >
{
  typedef DSG::Providers::Cube< GridType > GridProviderType;
  typedef Base< GridType > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef double            RangeFieldType;
  static const size_t       dimRange = 1;
  typedef DSG::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection > BoundaryInfoType;

  typedef Dune::Stuff::Functions::Constant
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ConstantFunctionType;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    ExpressionFunctionType;
  typedef Dune::Stuff::Functions::Spe10::Model1
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    Spe10Model1FunctionType;
  typedef Spe10Model1FunctionType DiffusionType;
  typedef ConstantFunctionType    ForceType;
  typedef ConstantFunctionType    DirichletType;
  typedef ConstantFunctionType    NeumannType;
  typedef ConstantFunctionType    ExactSolutionType;

  Spe10Model1(const size_t num_refinements = 1)
    : BaseType(create_initial_grid(), num_refinements)
    , boundary_info_()
    , diffusion_("perm_case1.dat",
                 {0.0, 0.0},
                 {5.0, 1.0})
    , force_(1)
    , dirichlet_(0)
    , neumann_(0)
  {}

  void print_header(std::ostream& out = std::cout) const
  {
    out << "+==========================================================+\n"
        << "|+========================================================+|\n"
        << "||  Testcase: SPE10, Model1                               ||\n"
        << "||  (see http://www.spe.org/web/csp/datasets/set01.htm)   ||\n"
        << "|+--------------------------------------------------------+|\n"
        << "||  domain = [0, 5] x [0, 1]                              ||\n"
        << "||  diffusion: spe10 model 1 scalar data                  ||\n"
        << "||  force     = 1                                         ||\n"
        << "||  dirichlet = 0                                         ||\n"
        << "||  reference solution: discrete solution on finest grid  ||\n"
        << "|+========================================================+|\n"
        << "+==========================================================+" << std::endl;
  } // ... print_header(...)

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const DiffusionType& diffusion() const
  {
    return diffusion_;
  }

  const ForceType& force() const
  {
    return force_;
  }

  const DirichletType& dirichlet() const
  {
    return dirichlet_;
  }

  const NeumannType& neumann() const
  {
    return neumann_;
  }

  bool provides_exact_solution() const
  {
    return false;
  }

  const ExactSolutionType& exact_solution() const
  {
    DUNE_THROW(Dune::InvalidStateException, "provides_exact_solution() == false");
    return neumann_;
  }

private:
  static std::unique_ptr< GridProviderType > create_initial_grid()
  {
    auto grid_provider = std::unique_ptr< GridProviderType >(new GridProviderType({0.0, 0.0},
                                                                                  {5.0, 1.0},
                                                                                  {100u, 20u}));
    grid_provider->grid().globalRefine(1);
    return grid_provider;
  } // ... create_initial_grid(...)

  const BoundaryInfoType boundary_info_;
  const DiffusionType diffusion_;
  const ForceType force_;
  const DirichletType dirichlet_;
  const NeumannType neumann_;
}; // class Spe10Model1


} // namespace EllipticTestCase


#if HAVE_ALUGRID

typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > AluConform2dGridType;

#define ALU_CONFORM_2D_ESTIMATOR_TESTCASES \
    EllipticTestCase::ESV07< AluConform2dGridType > \
  , EllipticTestCase::LocalThermalBlock< AluConform2dGridType > \
  , EllipticTestCase::Spe10Model1< AluConform2dGridType >

#define ALU_CONFORM_2D_TESTCASES \
    ALU_CONFORM_2D_ESTIMATOR_TESTCASES \
  , EllipticTestCase::ER07< AluConform2dGridType > \
  , EllipticTestCase::MixedBoundaryTypes< AluConform2dGridType > \

#endif // HAVE_ALUGRID


typedef testing::Types<
#if HAVE_ALUGRID
                      ALU_CONFORM_2D_ESTIMATOR_TESTCASES
#endif
                      > EllipticEstimatorTestCases;

typedef Dune::SGrid< 2, 2, double > S2dGridType;

typedef testing::Types<
#if HAVE_ALUGRID
                        ALU_CONFORM_2D_TESTCASES
#else
                        EllipticTestCase::LocalThermalBlock< S2dGridType >
                      , EllipticTestCase::Spe10Model1< S2dGridType >
                      , EllipticTestCase::ER07< S2dGridType >
                      , EllipticTestCase::MixedBoundaryTypes< S2dGridType >
#endif // HAVE_ALUGRID
                      > EllipticTestCases;


#endif // DUNE_GDT_TEST_ELLIPTIC_HH
