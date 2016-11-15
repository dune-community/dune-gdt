// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/projections.hh>
#include <dune/gdt/test/projections/base.hh>
#include <dune/gdt/test/grids.hh>

#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>
#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include <dune/gdt/test/projections.hh>
#include <dune/gdt/test/projections/l2.hh>

using namespace Dune::GDT::Test;

TEST_F(ProjectionTest, produces_correct_results)
{
  this->produces_correct_results();
}

TEST_F(L2ProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2ProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2ProjectionOperatorTest, free_function_callable)
{
  this->free_function_callable();
}
TEST_F(L2ProjectionOperatorTest, produces_correct_results)
{
  // RT : 0.0925927
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}

TEST_F(L2ProjectionLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2ProjectionLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2ProjectionLocalizableOperatorTest, produces_correct_results)
{
  // RT : 0.096226
  typedef typename TypeParam::GridViewType::Grid Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}
