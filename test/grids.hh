// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_GRIDS_HH
#define DUNE_GDT_TEST_GRIDS_HH

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/gdt/spaces/tools.hh>

using namespace Dune;
using namespace GDT;

#define SGRID_TYPES(dim) \
  typedef typename SpaceTools::LeafGridPartView<  SGrid< dim, dim >, false >::Type S ## dim ## dLeafGridPartType; \
  typedef typename SpaceTools::LevelGridPartView< SGrid< dim, dim >, false >::Type S ## dim ## dLevelGridPartType; \
  typedef typename SpaceTools::LeafGridPartView<  SGrid< dim, dim >, true  >::Type S ## dim ## dLeafGridViewType; \
  typedef typename SpaceTools::LevelGridPartView< SGrid< dim, dim >, true  >::Type S ## dim ## dLevelGridViewType;
SGRID_TYPES(1)
SGRID_TYPES(2)
SGRID_TYPES(3)
#undef SGRID_TYPES

#define YASPGRID_TYPES(dim) \
  typedef typename SpaceTools::LeafGridPartView<  YaspGrid< dim >, false >::Type Yasp ## dim ## dLeafGridPartType; \
  typedef typename SpaceTools::LevelGridPartView< YaspGrid< dim >, false >::Type Yasp ## dim ## dLevelGridPartType; \
  typedef typename SpaceTools::LeafGridPartView<  YaspGrid< dim >, true  >::Type Yasp ## dim ## dLeafGridViewType; \
  typedef typename SpaceTools::LevelGridPartView< YaspGrid< dim >, true  >::Type Yasp ## dim ## dLevelGridViewType;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES

#if HAVE_ALUGRID

typedef ALUGrid< 2, 2, simplex, conforming > AluConform2dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluConform2dGridType, false >::Type AluConform2dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluConform2dGridType, false >::Type AluConform2dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluConform2dGridType, true  >::Type AluConform2dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluConform2dGridType, true  >::Type AluConform2dLevelGridViewType;

typedef ALUGrid< 3, 3, simplex, conforming > AluConform3dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluConform3dGridType, false >::Type AluConform3dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluConform3dGridType, false >::Type AluConform3dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluConform3dGridType, true  >::Type AluConform3dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluConform3dGridType, true  >::Type AluConform3dLevelGridViewType;

typedef ALUGrid< 2, 2, simplex, nonconforming > AluSimplex2dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluSimplex2dGridType, false >::Type AluSimplex2dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluSimplex2dGridType, false >::Type AluSimplex2dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluSimplex2dGridType, true  >::Type AluSimplex2dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluSimplex2dGridType, true  >::Type AluSimplex2dLevelGridViewType;

typedef ALUGrid< 3, 3, simplex, nonconforming > AluSimplex3dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluSimplex3dGridType, false >::Type AluSimplex3dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluSimplex3dGridType, false >::Type AluSimplex3dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluSimplex3dGridType, true  >::Type AluSimplex3dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluSimplex3dGridType, true  >::Type AluSimplex3dLevelGridViewType;

typedef ALUGrid< 2, 2, cube, nonconforming > AluCube2dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluCube2dGridType, false >::Type AluCube2dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluCube2dGridType, false >::Type AluCube2dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluCube2dGridType, true  >::Type AluCube2dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluCube2dGridType, true  >::Type AluCube2dLevelGridViewType;

typedef ALUGrid< 3, 3, cube, nonconforming > AluCube3dGridType;
typedef typename SpaceTools::LeafGridPartView<  AluCube3dGridType, false >::Type AluCube3dLeafGridPartType;
typedef typename SpaceTools::LevelGridPartView< AluCube3dGridType, false >::Type AluCube3dLevelGridPartType;
typedef typename SpaceTools::LeafGridPartView<  AluCube3dGridType, true  >::Type AluCube3dLeafGridViewType;
typedef typename SpaceTools::LevelGridPartView< AluCube3dGridType, true  >::Type AluCube3dLevelGridViewType;

#endif // HAVE_ALUGRID

#endif // DUNE_GDT_TEST_GRIDS_HH
