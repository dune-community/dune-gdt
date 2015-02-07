// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_DARCY_HH
#define DUNE_GDT_TEST_OPERATORS_DARCY_HH

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/operators/darcy.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/fv/default.hh>

using namespace Dune;
using namespace GDT;


/**
 * \note This test assumes that DiscreteFunction, Operators::L2Projection, Products::L2, Products::H1Semi,
 *       Spaces::CG::FemBased, Spaces::RT::PdelabBased and Spaces::FV::Default work correctly.
 */
template< class SpaceTypes >
struct DarcyOperator
  : public ::testing::Test
{
  typedef typename SpaceTypes::first_type SourceSpaceType;
  typedef typename SpaceTypes::second_type RangeSpaceType;

  typedef typename RangeSpaceType::GridViewType     GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const size_t                   dimDomain = SourceSpaceType::dimDomain;
  typedef double RangeFieldType;

  typedef typename Dune::Stuff::LA::Container< RangeFieldType >::VectorType VectorType;

  void produces_correct_results() const
  {
    GridProviderType grid_provider(0.0, 1.0, 4);
    auto& grid = grid_provider.grid();
    grid.globalRefine(1);

    typedef Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > FunctionType;
    const FunctionType source("x", "x[0] * x[1]", 2, "source", {{"x[1]", "x[0]"}});

    const RangeSpaceType range_space(SpaceTools::GridPartView< RangeSpaceType >::create_leaf(grid));
    VectorType range_vector(range_space.mapper().size());
    DiscreteFunction< RangeSpaceType, VectorType > range(range_space, range_vector);

    const FunctionType function("x", "-1.0", 0);
    const Operators::Darcy< GridViewType, FunctionType > darcy_operator(range_space.grid_view(), function);
    darcy_operator.apply(source, range);

    const Stuff::Functions::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain >
      desired_output("x", std::vector< std::string >({"x[1]", "x[0]"}), 1,
                     "desired output",
                     {{"0.0", "1.0"}, {"1.0", "0.0"}});

    const Products::L2< GridViewType > l2_product(range_space.grid_view());
    const RangeFieldType l2_error = l2_product.induced_norm(desired_output - range);
    const RangeFieldType l2_error_expected = expected_result_("l2", desired_output, range_space.grid_view());
    EXPECT_LE(l2_error, l2_error_expected);

    const Products::H1Semi< GridViewType > h1_semi_product(range_space.grid_view());
    const RangeFieldType h1_error = h1_semi_product.induced_norm(desired_output - range);
    const RangeFieldType h1_error_expected = expected_result_("h1", desired_output, range_space.grid_view());
    EXPECT_LE(h1_error, h1_error_expected);
  } // ... produces_correct_results()

  template< class FunctionType, class GV >
  RangeFieldType expected_result_(const std::string type,
                                  const FunctionType& desired_output,
                                  const GV& grid_view) const
  {
    typedef typename SpaceTools::LeafGridPartView< GridType, RangeSpaceType::needs_grid_view >::Type GPV;
    if (std::is_base_of< Spaces::CG::FemBased< GPV, 1, RangeFieldType, dimDomain >
                       , RangeSpaceType >::value) {
      if (type == "l2")
        return 2.18e-16;
      else if (type == "h1")
        return 3.12e-15;
      else
        DUNE_THROW(Dune::Stuff::Exceptions::internal_error, type);
    } else if (std::is_base_of< Spaces::RT::PdelabBased< GPV, 0, RangeFieldType, dimDomain >
                              , RangeSpaceType >::value) {
      typedef Spaces::FV::Default< GV, RangeFieldType, dimDomain > FvSpaceType;
      const FvSpaceType fv_space(grid_view);
      VectorType fv_desired_output_vector(fv_space.mapper().size());
      DiscreteFunction< FvSpaceType, VectorType > fv_desired_output(fv_space, fv_desired_output_vector);
      const Operators::L2Projection< GV > l2_projection(grid_view);
      l2_projection.apply(desired_output, fv_desired_output);
      const Products::L2< GV > l2_product(grid_view);
      const Products::H1Semi< GV > h1_semi_product(grid_view);
      if (type == "l2")
        return 2.0 * l2_product.induced_norm(desired_output - fv_desired_output);
      else if (type == "h1")
        return h1_semi_product.induced_norm(desired_output - fv_desired_output);
      else
        DUNE_THROW(Dune::Stuff::Exceptions::internal_error, type);
    } else
      DUNE_THROW(Dune::Stuff::Exceptions::internal_error, type);
  } // ... expected_result_(...)
}; // struct DarcyOperator



#endif // DUNE_GDT_TEST_OPERATORS_DARCY_HH
