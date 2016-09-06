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
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
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

#if HAVE_ALUGRID

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> AluConform2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, false>::Type AluConform2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, true>::Type AluConform2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, true>::Type AluConform2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming> AluConform3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform3dGridType, false>::Type AluConform3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform3dGridType, false>::Type AluConform3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform3dGridType, true>::Type AluConform3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform3dGridType, true>::Type AluConform3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> AluSimplex2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming> AluSimplex3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> AluCube2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, false>::Type AluCube2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, false>::Type AluCube2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, true>::Type AluCube2dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, true>::Type AluCube2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming> AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, false>::Type AluCube3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, false>::Type AluCube3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, true>::Type AluCube3dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, true>::Type AluCube3dLevelGridViewType;

#endif // HAVE_ALUGRID

#endif // DUNE_GDT_TEST_GRIDS_HH
