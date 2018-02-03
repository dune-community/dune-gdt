// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/gdt/test/projections/lagrange.hh>
#include <dune/gdt/test/projections/l2-global.hh>
#include <dune/gdt/test/spaces/cg.hh>

using namespace Dune::GDT::Test;

// clang-format off
{% for Space, Name in config.spaces_with_names %}

typedef LagrangeProjectionOperatorTest<{{Space}}> LagrangeProjectionOperatorTest_{{Name}};
TEST_F(LagrangeProjectionOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(LagrangeProjectionOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(LagrangeProjectionOperatorTest_{{Name}}, free_function_callable)
{
  this->free_function_callable();
}
TEST_F(LagrangeProjectionOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results();
}

typedef LagrangeProjectionLocalizableOperatorTest<{{Space}}> LagrangeProjectionLocalizableOperatorTest_{{Name}};
TEST_F(LagrangeProjectionLocalizableOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(LagrangeProjectionLocalizableOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(LagrangeProjectionLocalizableOperatorTest_{{Name}}, produces_correct_results)
{
  this->produces_correct_results(this->default_tolerance);
}

typedef L2GlobalProjectionLocalizableOperatorTest<{{Space}}> L2GlobalProjectionLocalizableOperatorTest_{{Name}};
TEST_F(L2GlobalProjectionLocalizableOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2GlobalProjectionLocalizableOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2GlobalProjectionLocalizableOperatorTest_{{Name}}, produces_correct_results)
{

  using Grid = Dune::XT::Grid::extract_grid_t<typename L2GlobalProjectionLocalizableOperatorTest_{{Name}}::GridLayerType>;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}

typedef L2GlobalProjectionOperatorTest<{{Space}}> L2GlobalProjectionOperatorTest_{{Name}};
TEST_F(L2GlobalProjectionOperatorTest_{{Name}}, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TEST_F(L2GlobalProjectionOperatorTest_{{Name}}, constructible_by_factory)
{
  this->constructible_by_factory();
}
TEST_F(L2GlobalProjectionOperatorTest_{{Name}}, produces_correct_results)
{
  typedef Dune::XT::Grid::extract_grid_t<L2GlobalProjectionOperatorTest_{{Name}}::GridLayerType> Grid;
  const auto tolerance = Dune::XT::Grid::is_alugrid<Grid>::value ? this->alugrid_tolerance : this->default_tolerance;
  this->produces_correct_results(tolerance);
  this->produces_correct_results(tolerance);
}


{% endfor %}
// clang-format on
