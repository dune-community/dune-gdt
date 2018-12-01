// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_SPACES_RT_HH
#define DUNE_GDT_TEST_SPACES_RT_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/rt/default.hh>

#include <dune/gdt/test/grids.hh>

#include <dune/xt/grid/type_traits.hh>

#define SPACE_RT_YASPGRID(dd) Dune::GDT::RaviartThomasSpace<Yasp##dd##dLeafGridViewType, 0, double>

#define SPACES_RT SPACE_RT_YASPGRID(2), SPACE_RT_YASPGRID(3)

#define SPACE_RT_YASPGRID_LEVEL(dd) Dune::GDT::RaviartThomasSpace<Yasp##dd##dLevelGridViewType, 0, double>

#define SPACES_RT_LEVEL SPACE_RT_YASPGRID_LEVEL(2), SPACE_RT_YASPGRID_LEVEL(3)

#if HAVE_DUNE_ALUGRID

#  define SPACE_RT_ALUCONFORMGRID(dd) Dune::GDT::RaviartThomasSpace<AluConform##dd##dLeafGridViewType, 0, double>

#  define SPACE_RT_ALUCUBEGRID(dd) Dune::GDT::RaviartThomasSpace<AluCube##dd##dLeafGridViewType, 0, double>

#  define SPACES_RT_ALUGRID                                                                                            \
    SPACE_RT_ALUCONFORMGRID(2)                                                                                         \
    , SPACE_RT_ALUCUBEGRID(2), SPACE_RT_ALUCUBEGRID(3)


#  define SPACE_RT_ALUCONFORMGRID_LEVEL(dd) Dune::GDT::RaviartThomasSpace<AluConform##dd##dLevelGridViewType, 0, double>

#  define SPACE_RT_ALUCUBEGRID_LEVEL(dd) Dune::GDT::RaviartThomasSpace<AluCube##dd##dLevelGridViewType, 0, double>

#  define SPACES_RT_ALUGRID_LEVEL                                                                                      \
    SPACE_RT_ALUCONFORMGRID_LEVEL(2)                                                                                   \
    , SPACE_RT_ALUCUBEGRID_LEVEL(2), SPACE_RT_ALUCUBEGRID_LEVEL(3)

#endif // HAVE_DUNE_ALUGRID


#endif // DUNE_GDT_TEST_SPACES_RT_HH
