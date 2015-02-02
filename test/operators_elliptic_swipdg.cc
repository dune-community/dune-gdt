// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/main.hxx>

#if HAVE_DUNE_FEM && HAVE_EIGEN

#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/playground/operators/elliptic-swipdg.hh>

using namespace Dune;
using namespace Dune::GDT;


TEST(EllipticSWIPDGOperator, is_affinely_decomposable)
{
  static const unsigned int d = 2;
  typedef SGrid< d, d > GridType;
  typedef GridType::template Codim< 0 >::Entity E;
  typedef GridType::ctype D;
  typedef double R;
  static const unsigned int r = 1;
  auto grid_provider = Stuff::Grid::Providers::Cube< GridType >::create();
  typedef Spaces::DiscontinuousLagrangeProvider< GridType,
                                                 Stuff::Grid::ChooseLayer::leaf,
                                                 ChooseSpaceBackend::fem,
                                                 1, R, r > SpaceProvider;
  typedef SpaceProvider::Type SpaceType;
  auto space = SpaceProvider::create(*grid_provider);

  typedef typename SpaceType::GridViewType GridViewType;
  auto boundary_info = Stuff::Grid::BoundaryInfos::AllDirichlet< typename GridViewType::Intersection >::create();

  typedef Stuff::Functions::Constant< E, D, d, R, r >    ScalarFunctionType;
  typedef Stuff::Functions::Constant< E, D, d, R, d, d > TensorFunctionType;
  ScalarFunctionType one(17);
  TensorFunctionType tensor(42);

  auto two = Stuff::Functions::make_sum(one, one);

  typedef Stuff::LA::Container< R >::MatrixType MatrixType;

  auto one_op = Operators::make_elliptic_swipdg(one, tensor, *boundary_info, MatrixType(), space);
  auto two_op = Operators::make_elliptic_swipdg(*two, tensor, *boundary_info, MatrixType(), space);

  one_op->add(*two_op);
  one_op->assemble();

  auto tmp = one_op->matrix().copy();
  tmp.backend() += one_op->matrix().backend();

  tmp.backend() -= two_op->matrix().backend();
  EXPECT_EQ(0.0, tmp.sup_norm());
} // TEST(EllipticSWIPDGOperator, is_affinely_decomposable)

#else // HAVE_DUNE_FEM && HAVE_EIGEN

TEST(DISABLED_EllipticSWIPDGOperator, is_affinely_decomposable) {}

#endif
