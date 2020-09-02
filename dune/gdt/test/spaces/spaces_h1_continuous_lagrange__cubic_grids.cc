// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "h1_continuous_lagrange.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class G>
using Order2CubicContinuousLagrangeSpace = ContinuousLagrangeSpaceOnCubicLeafViewTest<G, 1, double, 2>;
TYPED_TEST_SUITE(Order2CubicContinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order2CubicContinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order2CubicContinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order3CubicContinuousLagrangeSpace = ContinuousLagrangeSpaceOnCubicLeafViewTest<G, 1, double, 3>;
TYPED_TEST_SUITE(Order3CubicContinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order3CubicContinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order3CubicContinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}
