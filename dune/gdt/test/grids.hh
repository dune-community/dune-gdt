// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_GRIDS_HH
#define DUNE_GDT_TEST_GRIDS_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/declaration.hh>
#endif

#include <dune/gdt/spaces/tools.hh>

#define YASPGRID_TYPES(dim)                                                                                            \
  typedef typename Dune::GDT::SpaceTools::                                                                             \
      LeafGridPartView<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>, false>::Type              \
          Yasp##dim##dLeafGridPartType;                                                                                \
  typedef typename Dune::GDT::SpaceTools::                                                                             \
      LevelGridPartView<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>, false>::Type             \
          Yasp##dim##dLevelGridPartType;                                                                               \
  typedef typename Dune::GDT::SpaceTools::                                                                             \
      LeafGridPartView<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>, true>::Type               \
          Yasp##dim##dLeafGridViewType;                                                                                \
  typedef typename Dune::GDT::SpaceTools::                                                                             \
      LevelGridPartView<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>, true>::Type              \
          Yasp##dim##dLevelGridViewType;                                                                               \
  typedef Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>> Yasp##dim##Grid;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES

#if HAVE_DUNE_ALUGRID
#if ALU3DGRID_PARALLEL
typedef Dune::ALUGridMPIComm AluComm;
#else
typedef Dune::ALUGridNoComm AluComm;
#endif

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming, AluComm> AluConform2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, false>::Type AluConform2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, true>::Type AluConform2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, true>::Type AluConform2dLevelGridViewType;

// typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming, AluComm> AluConform3dGridType;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>> AluConform3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform3dGridType, false>::Type AluConform3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform3dGridType, false>::Type AluConform3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform3dGridType, true>::Type AluConform3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform3dGridType, true>::Type AluConform3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming, AluComm> AluSimplex2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming, AluComm> AluSimplex3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming, AluComm> AluCube2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, false>::Type AluCube2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, false>::Type AluCube2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, true>::Type AluCube2dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, true>::Type AluCube2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, AluComm> AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, false>::Type AluCube3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, false>::Type AluCube3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, true>::Type AluCube3dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, true>::Type AluCube3dLevelGridViewType;

#endif // HAVE_DUNE_ALUGRID

template <class T>
double pdelab_rt_tolerance()
{
  typedef typename T::GridViewType::Grid Grid;
  constexpr auto dim = Grid::dimension;
  constexpr auto tolerance =
      Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.1 : 1.06) : (dim == 3 ? 2.05e-1 : 1.45e-1);
  return tolerance;
}

template <class T>
double pdelab_cg_tolerance()
{
  typedef typename T::GridViewType::Grid Grid;
  const auto dim = Grid::dimension;
  const auto tolerance = Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.35e-13 : 1.4e-14)
                                                                            : (dim == 3 ? 2.49e-14 : 1e-15);
  return tolerance;
}

template <class T>
double fem_cg_tolerance()
{
  typedef typename T::GridViewType::Grid Grid;
  const auto dim = Grid::dimension;
  const auto tolerance =
      Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.1e-13 : 1e-15) : (dim == 3 ? 2.49e-14 : 1e-15);
  return tolerance;
}
#endif // DUNE_GDT_TEST_GRIDS_HH
