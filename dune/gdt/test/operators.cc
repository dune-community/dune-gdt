// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/test_common.hh>

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/operator/products.hh>

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::Stuff::GridProviderCube<Dune::GridSelector::GridType> GridProviderType;
typedef typename GridProviderType::GridType GridType;
typedef Dune::grid::Part::Leaf::Const<GridType> GridPartType;
typedef typename GridPartType::template Codim<0>::EntityType EntityType;
typedef typename GridPartType::IntersectionType IntersectionType;
typedef Dune::Stuff::GridboundaryAllDirichlet<IntersectionType> BoundaryInfoType;
typedef typename GridPartType::ctype DomainFieldType;
static const unsigned int dimDomain = GridPartType::dimension;
typedef double RangeFieldType;
static const unsigned int dimRange = 1;
typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
    FunctionType;
typedef Dune::Stuff::LA::EigenDenseVector<RangeFieldType> VectorType;

typedef testing::
    Types<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<GridPartType, 1, RangeFieldType, dimRange>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<GridPartType, 2, RangeFieldType, dimRange>,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, 1, RangeFieldType, dimRange>,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, 2, RangeFieldType, dimRange>,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, 1, RangeFieldType, dimRange>,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, 2, RangeFieldType, dimRange>>
        SpaceTypes;

template <class SpaceType>
struct L2ProjectionOperator : public ::testing::Test
{
  void check() const
  {
    // prepare
    GridProviderType grid_provider;
    auto grid = grid_provider.grid();
    grid->globalRefine(4);
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::ProjectionOperator::L2<GridPartType> l2_projection_operator(*grid_part);
    l2_projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference<FunctionType, DiscreteFunctionType> difference(function, discrete_function);
    const Dune::GDT::ProductOperator::L2<GridPartType> l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > 1e-14)
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
  }
};

TYPED_TEST_CASE(L2ProjectionOperator, SpaceTypes);
TYPED_TEST(L2ProjectionOperator, compiles)
{
  this->check();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
