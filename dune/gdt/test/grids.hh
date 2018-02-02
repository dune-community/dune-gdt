// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_TEST_GRIDS_HH
#define DUNE_GDT_TEST_GRIDS_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/declaration.hh>
#endif

#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/type_traits.hh>

#define YASPGRID_TYPES(dim)                                                                                            \
  typedef typename Dune::XT::Grid::Layer<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>,         \
                                         Dune::XT::Grid::Layers::leaf,                                                 \
                                         Dune::XT::Grid::Backends::view>::type Yasp##dim##dLeafGridViewType;           \
  typedef typename Dune::XT::Grid::Layer<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>>,         \
                                         Dune::XT::Grid::Layers::level,                                                \
                                         Dune::XT::Grid::Backends::view>::type Yasp##dim##dLevelGridViewType;          \
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
typedef typename Dune::XT::Grid::Layer<AluConform2dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluConform2dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluConform2dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluConform2dLevelGridViewType;

// typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming, AluComm> AluConform3dGridType;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>> AluConform3dGridType;
typedef typename Dune::XT::Grid::Layer<AluConform3dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluConform3dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluConform3dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluConform3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming, AluComm> AluSimplex2dGridType;
typedef typename Dune::XT::Grid::Layer<AluSimplex2dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluSimplex2dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluSimplex2dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluSimplex2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming, AluComm> AluSimplex3dGridType;
typedef typename Dune::XT::Grid::Layer<AluSimplex3dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluSimplex3dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluSimplex3dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluSimplex3dLevelGridViewType;

typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming, AluComm> AluCube2dGridType;
typedef typename Dune::XT::Grid::Layer<AluCube2dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluCube2dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluCube2dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluCube2dLevelGridViewType;

typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, AluComm> AluCube3dGridType;
typedef typename Dune::XT::Grid::Layer<AluCube3dGridType,
                                       Dune::XT::Grid::Layers::leaf,
                                       Dune::XT::Grid::Backends::view>::type AluCube3dLeafGridViewType;
typedef typename Dune::XT::Grid::Layer<AluCube3dGridType,
                                       Dune::XT::Grid::Layers::level,
                                       Dune::XT::Grid::Backends::view>::type AluCube3dLevelGridViewType;

#endif // HAVE_DUNE_ALUGRID

template <class T>
double rt_tolerance()
{
  using Grid = Dune::XT::Grid::extract_grid_t<typename T::GridLayerType>;
  constexpr auto dim = Grid::dimension;
  constexpr auto tolerance =
      Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.1 : 1.06) : (dim == 3 ? 2.05e-1 : 1.45e-1);
  return tolerance;
}

template <class T>
double cg_tolerance()
{
  using Grid = Dune::XT::Grid::extract_grid_t<typename T::GridLayerType>;
  const auto dim = Grid::dimension;
  const auto tolerance = Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.35e-13 : 1.4e-14)
                                                                            : (dim == 3 ? 2.49e-14 : 1e-15);
  return tolerance;
}

#endif // DUNE_GDT_TEST_GRIDS_HH
