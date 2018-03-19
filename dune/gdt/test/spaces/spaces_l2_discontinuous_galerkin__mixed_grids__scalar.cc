// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "l2_discontinuous_galerkin.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class G>
using Order0ScalarMixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafViewTest<G, 1, 0>;
TYPED_TEST_CASE(Order0ScalarMixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order0ScalarMixedDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order1ScalarMixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafViewTest<G, 1, 1>;
TYPED_TEST_CASE(Order1ScalarMixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order1ScalarMixedDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order2ScalarMixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafViewTest<G, 1, 2>;
TYPED_TEST_CASE(Order2ScalarMixedDiscontinuousLagrangeSpace, MixedGrids);
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order2ScalarMixedDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}
