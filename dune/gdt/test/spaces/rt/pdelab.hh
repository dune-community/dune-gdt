// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_RT_PDELAB_HH
#define DUNE_GDT_TEST_SPACES_RT_PDELAB_HH

#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

#include <dune/gdt/test/grids.hh>


#define SPACE_RT_PDELAB_SGRID(dd) Dune::GDT::DunePdelabRtSpaceWrapper<S##dd##dLeafGridViewType, 0, double, dd>

#define SPACE_RT_PDELAB_YASPGRID(dd) Dune::GDT::DunePdelabRtSpaceWrapper<Yasp##dd##dLeafGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB                                                                                               \
  SPACE_RT_PDELAB_SGRID(2)                                                                                             \
  , SPACE_RT_PDELAB_SGRID(3), SPACE_RT_PDELAB_YASPGRID(2), SPACE_RT_PDELAB_YASPGRID(3)


#define SPACE_RT_PDELAB_SGRID_LEVEL(dd) Dune::GDT::DunePdelabRtSpaceWrapper<S##dd##dLevelGridViewType, 0, double, dd>

#define SPACE_RT_PDELAB_YASPGRID_LEVEL(dd)                                                                             \
  Dune::GDT::DunePdelabRtSpaceWrapper<Yasp##dd##dLevelGridViewType, 0, double, dd>

#define SPACES_RT_PDELAB_LEVEL                                                                                         \
  SPACE_RT_PDELAB_SGRID_LEVEL(2)                                                                                       \
  , SPACE_RT_PDELAB_SGRID_LEVEL(3), SPACE_RT_PDELAB_YASPGRID_LEVEL(2), SPACE_RT_PDELAB_YASPGRID_LEVEL(3)

#if HAVE_ALUGRID

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

#endif // HAVE_ALUGRID

#endif // DUNE_GDT_TEST_SPACES_RT_PDELAB_HH
