// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>

#include <dune/common/exceptions.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/fem/gridpart/levelgridpart.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/operator/prolongations.hh>

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::Stuff::LA::EigenDenseVector<double> VectorType;

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

// +-------------------------+
// |  * to test the products |
// +-------------------------+

template <class SpaceType, class ProductType>
struct LocalizableProduct
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef typename Dune::GDT::LocalizableProductInterface<typename ProductType::Traits> InterfaceType;

  static void fulfills_interface()
  {
    // static tests
    // * of the derived type
    typedef typename ProductType::Traits Traits;
    typedef typename ProductType::GridPartType D_GridPartType;
    typedef typename ProductType::RangeType D_RangeType;
    typedef typename ProductType::SourceType D_SourceType;
    typedef typename ProductType::FieldType D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridPartType I_GridPartType;
    typedef typename InterfaceType::RangeType I_RangeType;
    typedef typename InterfaceType::SourceType I_SourceType;
    typedef typename InterfaceType::FieldType I_FieldType;
    static_assert(std::is_same<ProductType, derived_type>::value, "");
    static_assert(std::is_same<D_GridPartType, I_GridPartType>::value, "");
    static_assert(std::is_same<D_RangeType, I_RangeType>::value, "");
    static_assert(std::is_same<D_SourceType, I_SourceType>::value, "");
    static_assert(std::is_same<D_FieldType, D_FieldType>::value, "");
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const FunctionType function("x", "1.0", 0);
    ProductType product(*grid_part, function, function);
    // dynamic tests
    // * of the derived type
    const D_GridPartType& d_gp = product.grid_part();
    if (&d_gp != &(*grid_part))
      DUNE_THROW(Dune::Exception, "");
    const D_RangeType& d_r = product.range();
    if (&d_r != &function)
      DUNE_THROW(Dune::Exception, "");
    const D_SourceType& d_s = product.source();
    if (&d_s != &function)
      DUNE_THROW(Dune::Exception, "");
    D_FieldType d_a = product.apply2();
    // * of the derived type as the interface
    InterfaceType& i_product   = static_cast<InterfaceType&>(product);
    const I_GridPartType& i_gp = i_product.grid_part();
    if (&i_gp != &d_gp)
      DUNE_THROW(Dune::Exception, "");
    const I_RangeType& i_r = i_product.range();
    if (&i_r != &d_r)
      DUNE_THROW(Dune::Exception, "");
    const I_SourceType& i_s = i_product.source();
    if (&i_s != &d_s)
      DUNE_THROW(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2();
    if (i_a != d_a)
      DUNE_THROW(Dune::Exception, "");
  }
}; // struct LocalizableProduct


template <class SpaceType, class ProductType, class VectorType>
struct AssemblableProduct
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef typename Dune::GDT::AssemblableProductInterface<typename ProductType::Traits> InterfaceType;

  static void fulfills_interface()
  {
    // static tests
    // * of the derived type
    typedef typename ProductType::Traits Traits;
    typedef typename ProductType::GridPartType D_GridPartType;
    typedef typename ProductType::RangeSpaceType D_RangeSpaceType;
    typedef typename ProductType::SourceSpaceType D_SourceSpaceType;
    typedef typename ProductType::MatrixType D_MatrixType;
    typedef typename ProductType::FieldType D_FieldType;
    // * of the derived type as the interface
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridPartType I_GridPartType;
    typedef typename InterfaceType::RangeSpaceType I_RangeSpaceType;
    typedef typename InterfaceType::SourceSpaceType I_SourceSpaceType;
    typedef typename InterfaceType::MatrixType I_MatrixType;
    typedef typename InterfaceType::FieldType I_FieldType;
    static_assert(std::is_same<ProductType, derived_type>::value, "");
    static_assert(std::is_same<D_GridPartType, I_GridPartType>::value, "");
    static_assert(std::is_same<D_RangeSpaceType, I_RangeSpaceType>::value, "");
    static_assert(std::is_same<D_SourceSpaceType, I_SourceSpaceType>::value, "");
    static_assert(std::is_same<D_MatrixType, D_MatrixType>::value, "");
    static_assert(std::is_same<D_FieldType, D_FieldType>::value, "");
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const SpaceType space(grid_part);
    const VectorType vector(space.mapper().size(), 1.0);
    ProductType product(*grid_part, space, space);
    // dynamic tests
    // * of the derived type
    const D_GridPartType& d_gp = product.grid_part();
    if (&d_gp != &(*grid_part))
      DUNE_THROW(Dune::Exception, "");
    const D_RangeSpaceType& d_r = product.range_space();
    if (&d_r != &space)
      DUNE_THROW(Dune::Exception, "");
    const D_SourceSpaceType& d_s = product.source_space();
    if (&d_s != &space)
      DUNE_THROW(Dune::Exception, "");
    D_MatrixType& d_m = product.matrix();
    D_FieldType d_a   = product.apply2(vector, vector);
    // * of the derived type as the interface
    InterfaceType& i_product   = static_cast<InterfaceType&>(product);
    const I_GridPartType& i_gp = i_product.grid_part();
    if (&i_gp != &d_gp)
      DUNE_THROW(Dune::Exception, "");
    const I_RangeSpaceType& i_r = i_product.range_space();
    if (&i_r != &d_r)
      DUNE_THROW(Dune::Exception, "");
    const I_SourceSpaceType& i_s = i_product.source_space();
    if (&i_s != &d_s)
      DUNE_THROW(Dune::Exception, "");
    I_MatrixType& i_m = i_product.matrix();
    if (&i_m != &d_m)
      DUNE_THROW(Dune::Exception, "");
    I_FieldType i_a = i_product.apply2(vector, vector);
    if (i_a != d_a)
      DUNE_THROW(Dune::Exception, "");
  }
}; // struct AssemblableProduct


template <class SpaceType>
struct GenericL2ProductOperator : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const Dune::GDT::ProductOperator::L2Generic<GridPartType> l2_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    auto l2_product      = l2_product_operator.apply2(function_1, function_1);
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    l2_product = l2_product_operator.apply2(function_2, function_2);
    error = l2_product - RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 3.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    l2_product = l2_product_operator.apply2(function_3, function_3);
    error = l2_product - RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 5.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
  }
}; // GenericL2ProductOperator

template <class SpaceType>
struct L2LocalizableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef Dune::GDT::ProductOperator::L2Localizable<GridPartType, FunctionType, FunctionType> ProductType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    ProductType l2_product_operator_1(*grid_part, function_1, function_1);
    auto l2_product      = l2_product_operator_1.apply2();
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    ProductType l2_product_operator_2(*grid_part, function_2, function_2);
    l2_product = l2_product_operator_2.apply2();
    error = l2_product - RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 3.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    ProductType l2_product_operator_3(*grid_part, function_3, function_3);
    l2_product = l2_product_operator_3.apply2();
    error = l2_product - RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0 / 5.0) << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    LocalizableProduct<SpaceType, ProductType>::fulfills_interface();
  }
}; // L2LocalizableProduct

template <class SpaceType>
struct L2AssemblableProduct : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<double> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<double> MatrixType;
  //  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  //  typedef typename SpaceType::DomainFieldType DomainFieldType;
  //  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  //  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  //  static const unsigned int                   dimRange = SpaceType::dimRange;

  typedef Dune::GDT::ProductOperator::L2Assemblable<GridPartType, SpaceType, SpaceType, MatrixType> ProductType;

  void produces_correct_results() const
  {
    //    // prepare
    //    const GridProviderType grid_provider(0.0, 1.0, 4u);
    //    const auto grid = grid_provider.grid();
    //    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    //    // test 1 (constant)
    //    const FunctionType function_1("x", "1.0", 0);
    //    ProductType l2_product_operator_1(*grid_part, function_1, function_1);
    //    auto l2_product = l2_product_operator_1.apply2();
    //    RangeFieldType error = l2_product - RangeFieldType(1.0);
    //    if (error > RangeFieldType(1e-15))
    //      DUNE_THROW(errors_are_not_as_expected,
    //                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
    //                 << " (difference: " << std::scientific << error << ")");
    //    // test 2 (linear)
    //    const FunctionType function_2("x", "x[0] - 1.0", 1);
    //    ProductType l2_product_operator_2(*grid_part, function_2, function_2);
    //    l2_product = l2_product_operator_2.apply2();
    //    error = l2_product - RangeFieldType(1.0/3.0);
    //    if (error > RangeFieldType(1e-15))
    //      DUNE_THROW(errors_are_not_as_expected,
    //                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
    //                 << " (difference: " << std::scientific << error << ")");
    //    // test 3 (quadratic)
    //    const FunctionType function_3("x", "x[0]*x[0]", 2);
    //    ProductType l2_product_operator_3(*grid_part, function_3, function_3);
    //    l2_product = l2_product_operator_3.apply2();
    //    error = l2_product - RangeFieldType(1.0/5.0);
    //    if (error > RangeFieldType(1e-15))
    //      DUNE_THROW(errors_are_not_as_expected,
    //                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
    //                 << " (difference: " << std::scientific << error << ")");
  } // ... produces_correct_results()

  void fulfills_interface() const
  {
    AssemblableProduct<SpaceType, ProductType, VectorType>::fulfills_interface();
  }
}; // L2AssemblableProduct

template <class SpaceType>
struct H1SemiProductOperator : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const Dune::GDT::ProductOperator::H1Semi<GridPartType> h1semi_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    auto h1semi_product  = h1semi_product_operator.apply2(function_1, function_1);
    RangeFieldType error = h1semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain * RangeFieldType(1.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 2 (linear)
    const FunctionType function_2(
        "x", "fake_value", 2, "affine gradient", {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    h1semi_product = h1semi_product_operator.apply2(function_2, function_2);
    error = h1semi_product - dimDomain * RangeFieldType(1.0 / 3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain * RangeFieldType(1.0 / 3.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
    // test 3 (quadratic)
    const FunctionType function_3(
        "x", "fake_value", 3, ", quadratic gradient", {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    h1semi_product = h1semi_product_operator.apply2(function_3, function_3);
    error = h1semi_product - dimDomain * RangeFieldType(1.0 / 5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain * RangeFieldType(1.0 / 5.0)
                                        << " (difference: "
                                        << std::scientific
                                        << error
                                        << ")");
  }
}; // H1SemiProductOperator

// +-------------------------------------+
// |  * to test the projection operators |
// +-------------------------------------+

template <class SpaceType, class ProjectionOperatorType>
struct ProjectionOperatorBase
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid      = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const ProjectionOperatorType projection_operator(*grid_part);
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference<FunctionType, DiscreteFunctionType> difference(function, discrete_function);
    const Dune::GDT::ProductOperator::L2Generic<GridPartType> l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // ProjectionOperatorType

template <class SpaceType>
struct L2ProjectionOperator
    : public ProjectionOperatorBase<SpaceType, Dune::GDT::ProjectionOperator::L2<typename SpaceType::GridPartType>>,
      public ::testing::Test
{
};

template <class SpaceType>
struct LagrangeProjectionOperator
    : public ProjectionOperatorBase<SpaceType,
                                    Dune::GDT::ProjectionOperator::Lagrange<typename SpaceType::GridPartType>>,
      public ::testing::Test
{
};

template <class SpaceType>
struct GenericProjectionOperator
    : public ProjectionOperatorBase<SpaceType,
                                    Dune::GDT::ProjectionOperator::Generic<typename SpaceType::GridPartType>>,
      public ::testing::Test
{
};

template <class SpaceType>
struct DirichletProjectionOperator : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 1u); // this has to be 1, otherwise the projection does not equal
    const auto grid      = grid_provider.grid(); // x[0] any more!
    const auto grid_part = std::make_shared<const GridPartType>(*grid);
    DomainType dirichlet_normal(0);
    dirichlet_normal[0] = DomainFieldType(1);
    const Dune::Stuff::GridboundaryNormalBased<typename GridPartType::IntersectionType> boundary_info(
        false, {dirichlet_normal});
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::ProjectionOperator::Dirichlet<GridPartType> projection_operator(*grid_part, boundary_info);
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference<FunctionType, DiscreteFunctionType> difference(function, discrete_function);
    const Dune::GDT::ProductOperator::L2Generic<GridPartType> l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // DirichletProjectionOperator

// +---------------------------------------+
// |  * to test the prolongation operators |
// +---------------------------------------+

template <class CoarseSpaceType, class FineSpaceType, class ProlongationOperatorType>
struct ProlongationOperatorBase
{
  typedef typename FineSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename FineSpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = FineSpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 2u);
    auto grid = grid_provider.grid();
    grid->globalRefine(1);
    const auto coarse_grid_part = std::make_shared<const GridPartType>(*grid, 0);
    assert(grid->maxLevel() > 0);
    const auto fine_grid_part = std::make_shared<const GridPartType>(*grid, grid->maxLevel());
    assert(fine_grid_part.size() > coarse_grid_part.size());
    // first, project an anlytical function onto the coarse grid
    const FunctionType function("x", "x[0]", 1, "function");
    const CoarseSpaceType coarse_space(coarse_grid_part);
    VectorType coarse_vector(coarse_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<CoarseSpaceType, VectorType> CoarseDiscreteFunctionType;
    CoarseDiscreteFunctionType coarse_discrete_function(coarse_space, coarse_vector, "coarse discrete function");
    const Dune::GDT::ProjectionOperator::Generic<GridPartType> coarse_projection_operator(*coarse_grid_part);
    coarse_projection_operator.apply(function, coarse_discrete_function);
    // since the projection operator was tested above we are confident this worked
    // but we check anyway (the L2 product operator was also tested above)
    const Dune::GDT::ProductOperator::L2Generic<GridPartType> coarse_l2_product_operator(*coarse_grid_part);
    const Dune::Stuff::Function::Difference<FunctionType, CoarseDiscreteFunctionType> coarse_difference(
        function, coarse_discrete_function);
    const auto coarse_l2_error = std::sqrt(coarse_l2_product_operator.apply2(coarse_difference, coarse_difference));
    if (coarse_l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "This should not happen, those operators were tested above!\n" << coarse_l2_error << " vs. "
                                                                                << RangeFieldType(1e-15));
    // now we prolong the discrete function from the coarse to the fine grid part
    const FineSpaceType fine_space(fine_grid_part);
    VectorType fine_vector(fine_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<FineSpaceType, VectorType> FineDiscreteFunctionType;
    FineDiscreteFunctionType fine_discrete_function(fine_space, fine_vector, "fine discrete function");
    const ProlongationOperatorType prolongation_operator(*fine_grid_part);
    prolongation_operator.apply(coarse_discrete_function, fine_discrete_function);
    // and measure the error
    const Dune::GDT::ProductOperator::L2Generic<GridPartType> fine_l2_product_operator(*fine_grid_part);
    const Dune::Stuff::Function::Difference<FunctionType, FineDiscreteFunctionType> fine_difference(
        function, fine_discrete_function);
    const auto fine_l2_error = std::sqrt(fine_l2_product_operator.apply2(fine_difference, fine_difference));
    if (fine_l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected, "\n" << fine_l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // ProlongationOperatorBase

template <class P>
struct L2ProlongationOperator
    : public ProlongationOperatorBase<typename P::first_type, typename P::second_type,
                                      Dune::GDT::ProlongationOperator::L2<typename P::second_type::GridPartType>>,
      public ::testing::Test
{
};

template <class P>
struct LagrangeProlongationOperator
    : public ProlongationOperatorBase<typename P::first_type, typename P::second_type,
                                      Dune::GDT::ProlongationOperator::Lagrange<typename P::second_type::GridPartType>>,
      public ::testing::Test
{
};

template <class P>
struct GenericProlongationOperator
    : public ProlongationOperatorBase<typename P::first_type, typename P::second_type,
                                      Dune::GDT::ProlongationOperator::Generic<typename P::second_type::GridPartType>>,
      public ::testing::Test
{
};

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

// +-----------------------------------------------------------------+
// | Therefore we need to define all grids and gridparts of interest |
// +-----------------------------------------------------------------+
#define SGRID_TYPES(dim)                                                                                               \
  typedef Dune::SGrid<dim, dim> S##dim##dGridType;                                                                     \
  typedef Dune::grid::Part::Leaf::Const<S##dim##dGridType> S##dim##dLeafGridPartType;                                  \
  typedef Dune::Fem::LevelGridPart<S##dim##dGridType> S##dim##dLevelGridPartType;
SGRID_TYPES(1)
SGRID_TYPES(2)
SGRID_TYPES(3)
#undef SGRID_TYPES

#define YASPGRID_TYPES(dim)                                                                                            \
  typedef Dune::YaspGrid<dim> Yasp##dim##dGridType;                                                                    \
  typedef Dune::grid::Part::Leaf::Const<Yasp##dim##dGridType> Yasp##dim##dLeafGridPartType;                            \
  typedef Dune::Fem::LevelGridPart<Yasp##dim##dGridType> Yasp##dim##dLevelGridPartType;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;
typedef Dune::grid::Part::Leaf::Const<AluConform2dGridType> AluConform2dLeafGridPartType;
typedef Dune::Fem::LevelGridPart<AluConform2dGridType> AluConform2dLevelGridPartType;
typedef Dune::ALUSimplexGrid<2, 2> AluSimplex2dGridType;
typedef Dune::grid::Part::Leaf::Const<AluSimplex2dGridType> AluSimplex2dLeafGridPartType;
typedef Dune::Fem::LevelGridPart<AluSimplex2dGridType> AluSimplex2dLevelGridPartType;
typedef Dune::ALUSimplexGrid<3, 3> AluSimplex3dGridType;
typedef Dune::grid::Part::Leaf::Const<AluSimplex3dGridType> AluSimplex3dLeafGridPartType;
typedef Dune::Fem::LevelGridPart<AluSimplex3dGridType> AluSimplex3dLevelGridPartType;
typedef Dune::ALUCubeGrid<3, 3> AluCube3dGridType;
typedef Dune::grid::Part::Leaf::Const<AluCube3dGridType> AluCube3dLeafGridPartType;
typedef Dune::Fem::LevelGridPart<AluCube3dGridType> AluCube3dLevelGridPartType;
#endif

// +----------------------------------------------------+
// |  * arguments for the product operator test structs |
// +----------------------------------------------------+

typedef testing::Types<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dLeafGridPartType, 1, double, 1>

                       ,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dLeafGridPartType, 1, double, 1>
#if HAVE_ALUGRID
                       ,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLeafGridPartType, 1, double, 1>,
                       Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dLeafGridPartType, 1, double, 1>
#endif
                       > ProductOperatorSpaceTypes;

// +-------------------------------------------------------+
// |  * arguments for the projection operator test structs |
// +-------------------------------------------------------+

#define L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID                                                                     \
  Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dLeafGridPartType, 1, double, 1>,         \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dLeafGridPartType, 1, double, 1>,     \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dLeafGridPartType, 1, double, 1>,     \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dLeafGridPartType, 2, double, 1>,     \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dLeafGridPartType, 2, double, 1>,     \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dLeafGridPartType, 2, double, 1>

typedef testing::Types<
#if HAVE_ALUGRID
    L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
    > L2ProjectionOperatorSpaceTypes;

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES                                                                       \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dLeafGridPartType, 1, double, 1>,                                   \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dLeafGridPartType, 1, double, 1>,                               \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dLeafGridPartType, 1, double, 1>                                \
                                                                                                                       \
      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dLeafGridPartType, 1, double, 1>,                          \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dLeafGridPartType, 1, double, 1>,                            \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dLeafGridPartType, 1, double, 1>

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID                                                               \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLeafGridPartType, 1, double, 1>,                          \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dLeafGridPartType, 1, double, 1>                          \
                                                                                                                       \
      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dLeafGridPartType, 1, double, 1>,      \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dLeafGridPartType, 1, double, 1>,        \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dLeafGridPartType, 1, double, 1>


typedef testing::Types<LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                       ,
                       LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                       > LagrangeProjectionOperatorSpaceTypes;

typedef testing::Types<LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                       ,
                       LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID, L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                       > GenericProjectionOperatorSpaceTypes;

typedef testing::Types<LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                       ,
                       LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                       > DirichletProjectionOperatorSpaceTypes;

#undef L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID

// +---------------------------------------------------------+
// |  * arguments for the prolongation operator test structs |
// +---------------------------------------------------------+

#define L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID                                                                   \
  /* all combinations which have DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper as FineSpaceType */              \
  std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>,               \
            Dune::GDT::DiscontinuousLagrangeSpace::                                                                    \
                FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,                                \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>                             \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::DiscontinuousLagrangeSpace::                                                              \
                      FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,                          \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>                             \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::DiscontinuousLagrangeSpace::                                                              \
                      FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,                          \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>

typedef testing::Types<
#if HAVE_ALUGRID
    L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
    > L2ProlongationOperatorSpaceTypes;

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES                                                                     \
  /* all combinations which have ContinuousLagrangeSpace::FemWrapper as FineSpaceType */                               \
  std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dLevelGridPartType, 1, double, 1>,                        \
            Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dLevelGridPartType, 1, double, 1>>,                       \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dLevelGridPartType, 1, double, 1>,                    \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dLevelGridPartType, 1, double, 1>>,                   \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dLevelGridPartType, 1, double, 1>,                    \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dLevelGridPartType, 1, double, 1>>                    \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dLevelGridPartType, 1, double, 1>,               \
                  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dLevelGridPartType, 1, double, 1>>,              \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dLevelGridPartType, 1, double, 1>,                 \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dLevelGridPartType, 1, double, 1>>,                \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dLevelGridPartType, 1, double, 1>,                 \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dLevelGridPartType, 1, double, 1>>

#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID                                                             \
  /* all combinations which have ContinuousLagrangeSpace::FemWrapper as FineSpaceType */                               \
  std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>,               \
            Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,              \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,          \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>>           \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,        \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,          \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>           \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,        \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,          \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>           \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dLevelGridPartType, 1, double, 1>,            \
                  Dune::GDT::ContinuousLagrangeSpace::                                                                 \
                      FemWrapper<AluCube3dLevelGridPartType,                                                           \
                                 1,                                                                                    \
                                 double,                                                                               \
                                 1>> /* all combinations which have                                                    \
                                         ContinuousLagrangeSpace::FemLocalfunctionsWrapper as                          \
                                         FineSpaceType */                                                              \
      ,                                                                                                                \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dLevelGridPartType, 1, double, 1>,           \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluConform2dLevelGridPartType, 1, double, 1>>                             \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::ContinuousLagrangeSpace::                                                                 \
                      FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,                          \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex2dLevelGridPartType, 1, double, 1>>                             \
                                                                                                                       \
      , std::pair<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,         \
                  Dune::GDT::ContinuousLagrangeSpace::                                                                 \
                      FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,                          \
      std::pair<Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>,                            \
      std::pair<Dune::GDT::DiscontinuousLagrangeSpace::                                                                \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>,                             \
                Dune::GDT::ContinuousLagrangeSpace::                                                                   \
                    FemLocalfunctionsWrapper<AluSimplex3dLevelGridPartType, 1, double, 1>>

typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                       ,
                       LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                       > LagrangeProlongationOperatorSpaceTypes;

typedef testing::Types<LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                       ,
                       LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID, L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                       > GenericProlongationOperatorSpaceTypes;

#undef L2_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES

// +--------------------------------------------------------------------------------------+
// | 3rd we combine all test structs with their appropriate arguments to create the tests |
// | (comment out the following lines if you do not want a test to be run)                |
// +--------------------------------------------------------------------------------------+

// +-------------------------------+
// |  * the product operator tests |
// +-------------------------------+

TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2LocalizableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2LocalizableProduct, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(L2AssemblableProduct, ProductOperatorSpaceTypes);
TYPED_TEST(L2AssemblableProduct, produces_correct_results)
{
  this->produces_correct_results();
}

TYPED_TEST_CASE(GenericL2ProductOperator, ProductOperatorSpaceTypes);
TYPED_TEST(GenericL2ProductOperator, produces_correct_results)
{
  this->produces_correct_results();
}

// TYPED_TEST_CASE(H1SemiProductOperator, ProductOperatorSpaceTypes);
// TYPED_TEST(H1SemiProductOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// +----------------------------------+
// |  * the projection operator tests |
// +----------------------------------+

// TYPED_TEST_CASE(L2ProjectionOperator, L2ProjectionOperatorSpaceTypes);
// TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(LagrangeProjectionOperator, LagrangeProjectionOperatorSpaceTypes);
// TYPED_TEST(LagrangeProjectionOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(GenericProjectionOperator, GenericProjectionOperatorSpaceTypes);
// TYPED_TEST(GenericProjectionOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(DirichletProjectionOperator, DirichletProjectionOperatorSpaceTypes);
// TYPED_TEST(DirichletProjectionOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// +------------------------------------+
// |  * the prolongation operator tests |
// +------------------------------------+

// TYPED_TEST_CASE(L2ProlongationOperator, L2ProlongationOperatorSpaceTypes);
// TYPED_TEST(L2ProlongationOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(LagrangeProlongationOperator, LagrangeProlongationOperatorSpaceTypes);
// TYPED_TEST(LagrangeProlongationOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// TYPED_TEST_CASE(GenericProlongationOperator, GenericProlongationOperatorSpaceTypes);
// TYPED_TEST(GenericProlongationOperator, produces_correct_results) {
//  this->produces_correct_results();
//}

// +--------------------------------------------------------------------------------------+
// | 4th we run all the tests                                                             |
// | (run the resulting executable with '--gtest_catch_exceptions=0' to see an exception) |
// +--------------------------------------------------------------------------------------+

int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
