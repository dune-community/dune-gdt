// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   Rene Milk       (2014, 2016, 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_OPERATORS_DARCY_HH
#define DUNE_GDT_TEST_OPERATORS_DARCY_HH

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/la/container.hh>
#include <dune/xt/functions/expression.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/operators/darcy.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/operators/laplace.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt/default.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \note This test assumes that DiscreteFunction, Operators::L2Projection, Products::L2, Products::H1Semi,
 *       ContinuousLagrangeSpace, RaviartThomasSpace and FvSpace work correctly.
 * \todo This test is rather old and could be refactored in terms of the other operator tests.
 * \todo Missing ctor and make_darcy_operator tests.
 */
template <class SpaceTypes>
struct DarcyOperatorTest : public ::testing::Test
{
  typedef typename SpaceTypes::first_type SourceSpaceType;
  typedef typename SpaceTypes::second_type RangeSpaceType;

  typedef typename RangeSpaceType::GridLayerType GridLayerType;
  typedef XT::Grid::extract_grid_t<GridLayerType> GridType;
  typedef XT::Grid::GridProvider<GridType> GridProviderType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = SourceSpaceType::dimDomain;
  typedef double RangeFieldType;

  typedef typename Dune::XT::LA::Container<RangeFieldType>::VectorType VectorType;

  void produces_correct_results() const
  {
    GridProviderType grid_provider(XT::Grid::make_cube_grid<GridType>(0.0, 1.0, 4));
    grid_provider.global_refine(1);

    typedef XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1> FunctionType;
    const FunctionType source("x", "x[0] * x[1]", 2, "source", {"x[1]", "x[0]"});

    const RangeSpaceType range_space(
        grid_provider.template layer<XT::Grid::Layers::leaf, RangeSpaceType::layer_backend>());
    VectorType range_vector(range_space.mapper().size());
    DiscreteFunction<RangeSpaceType, VectorType> range(range_space, range_vector);

    const FunctionType function("x", "-1.0", 0);
    const DarcyOperator<GridLayerType, FunctionType> darcy_operator(range_space.grid_layer(), function);
    darcy_operator.apply(source, range);

    const XT::Functions::ExpressionFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain>
        desired_output(
            "x", std::vector<std::string>({"x[1]", "x[0]"}), 1, "desired output", {{"0.0", "1.0"}, {"1.0", "0.0"}});

    const RangeFieldType l2_error = make_l2_operator(range_space.grid_layer(), 2)->induced_norm(desired_output - range);
    const RangeFieldType l2_error_expected = expected_result_("l2", desired_output, range_space.grid_layer());
    EXPECT_LE(l2_error, l2_error_expected);

    const RangeFieldType h1_error =
        make_laplace_operator(range_space.grid_layer(), 2)->induced_norm(desired_output - range);
    const RangeFieldType h1_error_expected = expected_result_("h1", desired_output, range_space.grid_layer());
    EXPECT_LE(h1_error, h1_error_expected);
  } // ... produces_correct_results()

  template <class FunctionType, class GL>
  RangeFieldType
  expected_result_(const std::string type, const FunctionType& desired_output, const GL& grid_layer) const
  {
    if (std::is_base_of<Dune::GDT::ContinuousLagrangeSpace<GL, 1, double>, RangeSpaceType>::value) {
      if (type == "l2")
        return 2.18e-16;
      else if (type == "h1")
        return 3.12e-15;
      else
        DUNE_THROW(Dune::XT::Common::Exceptions::internal_error, type);
    } else if (std::is_base_of<RaviartThomasSpace<GL, 0, RangeFieldType>, RangeSpaceType>::value) {
      typedef FvSpace<GL, RangeFieldType, dimDomain> FvSpaceType;
      const FvSpaceType fv_space(grid_layer);
      VectorType fv_desired_output_vector(fv_space.mapper().size());
      DiscreteFunction<FvSpaceType, VectorType> fv_desired_output(fv_space, fv_desired_output_vector);
      project(desired_output, fv_desired_output);
      if (type == "l2") {
        return 2.0 * make_l2_operator(grid_layer, 2)->induced_norm(desired_output - fv_desired_output);
      } else if (type == "h1") {
        return make_laplace_operator(grid_layer, 2)->induced_norm(desired_output - fv_desired_output);
      } else
        DUNE_THROW(Dune::XT::Common::Exceptions::internal_error, type);
    } else
      DUNE_THROW(Dune::XT::Common::Exceptions::internal_error, type);
  } // ... expected_result_(...)
}; // struct DarcyOperatorTest


} // namespace Test
} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_TEST_OPERATORS_DARCY_HH
