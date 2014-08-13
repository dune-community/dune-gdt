// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATOR_PRODUCTS
#define DUNE_GDT_TEST_OPERATOR_PRODUCTS

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>

#if HAVE_ALUGRID
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/grid/alugrid.hh>
#   include <dune/grid/sgrid.hh>
#   include <dune/grid/yaspgrid.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_ALUGRID

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/continuouslagrange/fem.hh>
#include <dune/gdt/playground/spaces/discontinuouslagrange/fem.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/operators/prolongations.hh>

typedef Dune::Stuff::LA::IstlDenseVector< double > VectorType;

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

// +-------------------------+
// |  * to test the products |
// +-------------------------+

//      - first the interfaces

template< class SpaceType, class ProductType >
struct ProductBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef typename Dune::GDT::ProductInterface< typename ProductType::Traits > InterfaceType;

  static void fulfills_interface()
  {
    // static tests
    // * of the derived type
#include <dune/stuff/common/disable_warnings.hh>
    typedef typename ProductType::Traits        Traits;
#include <dune/stuff/common/reenable_warnings.hh>
    typedef typename ProductType::GridViewType  D_GridViewType;
    typedef typename ProductType::FieldType     D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type  derived_type;
    typedef typename InterfaceType::GridViewType  I_GridViewType;
    typedef typename InterfaceType::FieldType     I_FieldType;
    static_assert(std::is_same< ProductType, derived_type >::value, "");
    static_assert(std::is_same< D_GridViewType, I_GridViewType >::value, "");
    static_assert(std::is_same< D_FieldType, D_FieldType >::value, "");
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "1.0", 0, "function", {{"1.0", "1.0", "1.0"}});
    ProductType product(*(space.grid_view()));
    // dynamic tests
    // * of the derived type
    const D_GridViewType& d_gp = product.grid_view();
    if (&d_gp != &(*(space.grid_view()))) DUNE_THROW(Dune::Exception, "");
    D_FieldType d_a = product.apply2(function, function);
    // * of the derived type as the interface
    InterfaceType& i_product = static_cast< InterfaceType& >(product);
    const I_GridViewType& i_gp = i_product.grid_view();
    if (&i_gp != &d_gp) DUNE_THROW(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2(function, function);
    if (Dune::Stuff::Common::FloatCmp::ne(i_a, d_a)) DUNE_THROW(Dune::Exception, "");
  }
}; // struct ProductBase

// +----------------------------------------------------+
// |  * arguments for the product operator test structs |
// +----------------------------------------------------+


#define SGRID_TYPES(dim) \
  typedef Dune::SGrid< dim, dim >                                 S ## dim ## dGridType; \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S ## dim ## dGridType, false >::Type  S ## dim ## dLeafGridPartType; \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< S ## dim ## dGridType, false >::Type S ## dim ## dLevelGridPartType;
SGRID_TYPES(1)
SGRID_TYPES(2)
SGRID_TYPES(3)
#undef SGRID_TYPES

#define YASPGRID_TYPES(dim) \
  typedef Dune::YaspGrid< dim >                                     Yasp ## dim ## dGridType; \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp ## dim ## dGridType, false >::Type   Yasp ## dim ## dLeafGridPartType; \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView< Yasp ## dim ## dGridType, false >::Type  Yasp ## dim ## dLevelGridPartType;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES

#if HAVE_ALUGRID
typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > AluConform2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, false >::Type   AluConform2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluConform2dGridType, false >::Type  AluConform2dLevelGridPartType;
typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming > AluSimplex2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex2dGridType, false >::Type   AluSimplex2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluSimplex2dGridType, false >::Type  AluSimplex2dLevelGridPartType;
typedef Dune::ALUGrid< 3, 3, Dune::simplex, Dune::nonconforming > AluSimplex3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex3dGridType, false >::Type   AluSimplex3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluSimplex3dGridType, false >::Type  AluSimplex3dLevelGridPartType;
typedef Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming > AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluCube3dGridType, false >::Type  AluCube3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView< AluCube3dGridType, false >::Type AluCube3dLevelGridPartType;
#endif

typedef testing::Types< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLeafGridPartType, 1, double, 1 >

                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLeafGridPartType, 1, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 >
                      , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLeafGridPartType, 1, double, 1 >
#endif
                      > ProductOperatorSpaceTypes;

#endif // DUNE_GDT_TEST_OPERATOR_PRODUCTS
