// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1

// this one has to come first
#include <dune/stuff/test/main.hxx>

#include "spaces_fv.hh"
#include "spaces_dg_fem.hh"
#include "spaces_cg_pdelab.hh"

#include "products_h1.hh"


typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1)> ConstantSpaces;

TYPED_TEST_CASE(H1SemiLocalizableProduct, ConstantSpaces);
TYPED_TEST(H1SemiLocalizableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(H1SemiLocalizableProduct, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(H1SemiLocalizableProduct, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(H1SemiLocalizableProduct, quadratic_arguments)
{
  this->quadratic_arguments();
}


#if HAVE_DUNE_FEM

typedef testing::Types<SPACE_DG_FEM_SGRID(1, 1, 3), SPACE_DG_FEM_SGRID(2, 1, 3), SPACE_DG_FEM_SGRID(3, 1, 3)>
    CubicSpaces;

TYPED_TEST_CASE(H1SemiAssemblableProduct, CubicSpaces);
TYPED_TEST(H1SemiAssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(H1SemiAssemblableProduct, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(H1SemiAssemblableProduct, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(H1SemiAssemblableProduct, quadratic_arguments)
{
  this->quadratic_arguments();
}

#elif HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

typedef testing::Types<SPACE_CG_PDELAB_SGRID(1, 1, 1), SPACE_CG_PDELAB_SGRID(2, 1, 1), SPACE_CG_PDELAB_SGRID(3, 1, 1)>
    LinearSpaces;

TYPED_TEST_CASE(H1SemiAssemblableProduct, LinearSpaces);
TYPED_TEST(H1SemiAssemblableProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(H1SemiAssemblableProduct, constant_arguments)
{
  this->constant_arguments();
}
TEST(DISABLED_H1SemiAssemblableProduct, linear_arguments)
{
}
TEST(DISABLED_H1SemiAssemblableProduct, quadratic_arguments)
{
}

#else // HAVE_DUNE_PDELAB // HAVE_DUNE_FEM

TEST(DISABLED_H1SemiAssemblableProduct, fulfills_interface)
{
}
TEST(DISABLED_H1SemiAssemblableProduct, constant_arguments)
{
}
TEST(DISABLED_H1SemiAssemblableProduct, linear_arguments)
{
}
TEST(DISABLED_H1SemiAssemblableProduct, quadratic_arguments)
{
}

#endif // HAVE_DUNE_PDELAB // HAVE_DUNE_FEM


TYPED_TEST_CASE(H1SemiProduct, ConstantSpaces);
TYPED_TEST(H1SemiProduct, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(H1SemiProduct, constant_arguments)
{
  this->constant_arguments();
}
TYPED_TEST(H1SemiProduct, linear_arguments)
{
  this->linear_arguments();
}
TYPED_TEST(H1SemiProduct, quadratic_arguments)
{
  this->quadratic_arguments();
}
