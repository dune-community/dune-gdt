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
using Order0ScalarSimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest<G, 1, 0>;
TYPED_TEST_CASE(Order0ScalarSimplicialDiscontinuousLagrangeSpace, SimplicialGridsForSpaceTest);
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order0ScalarSimplicialDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order1ScalarSimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest<G, 1, 1>;
TYPED_TEST_CASE(Order1ScalarSimplicialDiscontinuousLagrangeSpace, SimplicialGridsForSpaceTest);
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order1ScalarSimplicialDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order2ScalarSimplicialDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnSimplicialLeafViewTest<G, 1, 2>;
TYPED_TEST_CASE(Order2ScalarSimplicialDiscontinuousLagrangeSpace, SimplicialGridsForSpaceTest);
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order2ScalarSimplicialDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order0ScalarCubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, 1, 0>;
TYPED_TEST_CASE(Order0ScalarCubicDiscontinuousLagrangeSpace, CubicGridsForSpaceTest);
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order0ScalarCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order1ScalarCubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, 1, 1>;
TYPED_TEST_CASE(Order1ScalarCubicDiscontinuousLagrangeSpace, CubicGridsForSpaceTest);
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order1ScalarCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order2ScalarCubicDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnCubicLeafViewTest<G, 1, 2>;
TYPED_TEST_CASE(Order2ScalarCubicDiscontinuousLagrangeSpace, CubicGridsForSpaceTest);
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order2ScalarCubicDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order0ScalarPrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafViewTest<G, 1, 0>;
TYPED_TEST_CASE(Order0ScalarPrismDiscontinuousLagrangeSpace, PrismGridsForSpaceTest);
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order0ScalarPrismDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order1ScalarPrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafViewTest<G, 1, 1>;
TYPED_TEST_CASE(Order1ScalarPrismDiscontinuousLagrangeSpace, PrismGridsForSpaceTest);
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order1ScalarPrismDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order2ScalarPrismDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnPrismLeafViewTest<G, 1, 2>;
TYPED_TEST_CASE(Order2ScalarPrismDiscontinuousLagrangeSpace, PrismGridsForSpaceTest);
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, gives_correct_identification)
{
  this->gives_correct_identification();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_of_lagrange_space_exists_on_each_element_with_correct_order();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_of_lagrange_space_on_each_element();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, mapper_maps_correctly)
{
  this->mapper_of_discontinuous_space_maps_correctly();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, lagrange_points_exist_on_each_element_with_correct_size)
{
  this->lagrange_points_exist_on_each_element_with_correct_size();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, basis_is_lagrange_basis)
{
  this->basis_is_lagrange_basis();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, basis_jacobians_seem_to_be_correct)
{
  this->template basis_jacobians_of_lagrange_space_seem_to_be_correct<>();
}
TYPED_TEST(Order2ScalarPrismDiscontinuousLagrangeSpace, local_interpolation_seems_to_be_correct)
{
  this->local_interpolation_seems_to_be_correct();
}


template <class G>
using Order0ScalarMixedDiscontinuousLagrangeSpace = DiscontinuousLagrangeSpaceOnMixedLeafViewTest<G, 1, 0>;
TYPED_TEST_CASE(Order0ScalarMixedDiscontinuousLagrangeSpace, MixedGridsForSpaceTest);
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
TYPED_TEST_CASE(Order1ScalarMixedDiscontinuousLagrangeSpace, MixedGridsForSpaceTest);
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
TYPED_TEST_CASE(Order2ScalarMixedDiscontinuousLagrangeSpace, MixedGridsForSpaceTest);
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
