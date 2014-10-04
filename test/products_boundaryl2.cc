// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv_default.hh"
#include "spaces_dg_fem.hh"

#include "products_boundaryl2.hh"

typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(BoundaryL2LocalizableProduct, ConstantSpaces);
TYPED_TEST(BoundaryL2LocalizableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(BoundaryL2LocalizableProduct, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(BoundaryL2LocalizableProduct, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(BoundaryL2LocalizableProduct, quadratic_arguments)
{
  this->quadratic_arguments();
}


#if HAVE_DUNE_FEM

typedef testing::Types<SPACE_DG_FEM_SGRID(1, 1, 2), SPACE_DG_FEM_SGRID(2, 1, 2), SPACE_DG_FEM_SGRID(3, 1, 2)>
    QuadraticSpaces;

TYPED_TEST_CASE(BoundaryL2AssemblableProduct, QuadraticSpaces);
TYPED_TEST(BoundaryL2AssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(BoundaryL2AssemblableProduct, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(BoundaryL2AssemblableProduct, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(BoundaryL2AssemblableProduct, quadratic_arguments)
{
  this->quadratic_arguments();
}

#else // HAVE_DUNE_FEM

TYPED_TEST_CASE(BoundaryL2AssemblableProduct, ConstantSpaces);
TYPED_TEST(BoundaryL2AssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(BoundaryL2AssemblableProduct, constant_arguments)
{
  this->constant_arguments();
}
TEST(DISABLED_BoundaryL2AssemblableProduct, linear_arguments)
{
}
TEST(DISABLED_BoundaryL2AssemblableProduct, quadratic_arguments)
{
}

#endif // HAVE_DUNE_FEM


TYPED_TEST_CASE(BoundaryL2Product, ConstantSpaces);
TYPED_TEST(BoundaryL2Product, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(BoundaryL2Product, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(BoundaryL2Product, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(BoundaryL2Product, quadratic_arguments)
{
  this->quadratic_arguments();
}
