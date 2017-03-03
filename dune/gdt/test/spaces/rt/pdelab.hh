// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_TEST_SPACES_RT_PDELAB_HH
#define DUNE_GDT_TEST_SPACES_RT_PDELAB_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

#include <dune/gdt/test/grids.hh>

#include <dune/xt/grid/type_traits.hh>

#define SPACE_RT_PDELAB_YASPGRID(dd) Dune::GDT::DunePdelabRtSpaceWrapper<Yasp##dd##dLeafGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB SPACE_RT_PDELAB_YASPGRID(2), SPACE_RT_PDELAB_YASPGRID(3)

#define SPACE_RT_PDELAB_YASPGRID_LEVEL(dd)                                                                             \
  Dune::GDT::DunePdelabRtSpaceWrapper<Yasp##dd##dLevelGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB_LEVEL SPACE_RT_PDELAB_YASPGRID_LEVEL(2), SPACE_RT_PDELAB_YASPGRID_LEVEL(3)

#if HAVE_DUNE_ALUGRID

#define SPACE_RT_PDELAB_ALUCONFORMGRID(dd)                                                                             \
  Dune::GDT::DunePdelabRtSpaceWrapper<AluConform##dd##dLeafGridViewType, 0, double, dd>

#define SPACE_RT_PDELAB_ALUCUBEGRID(dd)                                                                                \
  Dune::GDT::DunePdelabRtSpaceWrapper<AluCube##dd##dLeafGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB_ALUGRID                                                                                       \
  SPACE_RT_PDELAB_ALUCONFORMGRID(2)                                                                                    \
  , SPACE_RT_PDELAB_ALUCUBEGRID(2), SPACE_RT_PDELAB_ALUCUBEGRID(3)


#define SPACE_RT_PDELAB_ALUCONFORMGRID_LEVEL(dd)                                                                       \
  Dune::GDT::DunePdelabRtSpaceWrapper<AluConform##dd##dLevelGridViewType, 0, double, dd>

#define SPACE_RT_PDELAB_ALUCUBEGRID_LEVEL(dd)                                                                          \
  Dune::GDT::DunePdelabRtSpaceWrapper<AluCube##dd##dLevelGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB_ALUGRID_LEVEL                                                                                 \
  SPACE_RT_PDELAB_ALUCONFORMGRID_LEVEL(2)                                                                              \
  , SPACE_RT_PDELAB_ALUCUBEGRID_LEVEL(2), SPACE_RT_PDELAB_ALUCUBEGRID_LEVEL(3)

#endif // HAVE_DUNE_ALUGRID


#endif // DUNE_GDT_TEST_SPACES_RT_PDELAB_HH
