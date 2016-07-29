// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/stuff/test/main.hxx>

#include "prolongations/l2.hh"
#include "spaces/cg/pdelab.hh"

using namespace Dune::GDT::Test;

#if HAVE_DUNE_PDELAB


typedef testing::Types<SPACES_CG_PDELAB_LEVEL(1)
#if HAVE_DUNE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID_LEVEL(1)
#endif
                       > SpaceTypes;

template <class T>
double get_tolerance(const T& param)
{
  typedef typename T::GridViewType::Grid Grid;
  const auto dim = param.dimDomain;
  const auto tolerance = Dune::Stuff::Grid::is_conforming_alugrid<Grid>::value
      ? (dim == 3 ? 2.49e-14 : 1e-15)
      : (dim == 3 ? 2.49e-14 : 1e-15);
  return tolerance;
}

TYPED_TEST_CASE(L2ProlongationLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2ProlongationLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2ProlongationLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2ProlongationLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results(get_tolerance(*this));
}


#else // HAVE_DUNE_PDELAB


TEST(DISABLED_L2ProlongationLocalizableOperatorTest, constructible_by_ctor)
{
}
TEST(DISABLED_L2ProlongationLocalizableOperatorTest, constructible_by_factory)
{
}
TEST(DISABLED_L2ProlongationLocalizableOperatorTest, produces_correct_results)
{
}


#endif // HAVE_DUNE_PDELAB
