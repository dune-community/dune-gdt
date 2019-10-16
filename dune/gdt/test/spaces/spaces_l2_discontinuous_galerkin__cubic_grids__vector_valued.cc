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

#include "l2_discontinuous_galerkin.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class G>
using Order0VectorValuedCubicDiscontinuousLagrangeSpace =
    DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, G::dimension + 1, double, 0>;
TYPED_TEST_CASE(Order0VectorValuedCubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
//  TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
//  {
//    this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
//  }
TYPED_TEST(Order0VectorValuedCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order1VectorValuedCubicDiscontinuousLagrangeSpace =
    DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, G::dimension + 1, double, 1>;
TYPED_TEST_CASE(Order1VectorValuedCubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
//  TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
//  {
//    this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
//  }
TYPED_TEST(Order1VectorValuedCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order2VectorValuedCubicDiscontinuousLagrangeSpace =
    DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, G::dimension + 1, double, 2>;
TYPED_TEST_CASE(Order2VectorValuedCubicDiscontinuousLagrangeSpace, CubicGrids);
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
//  TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
//  {
//    this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
//  }
TYPED_TEST(Order2VectorValuedCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}
