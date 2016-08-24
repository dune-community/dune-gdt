// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

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

template <class T>
double
pdelab_rt_tolerance(const T& param)
{
  typedef typename T::GridViewType::Grid Grid;
  const auto dim = param.dimDomain;
  const auto tolerance =
      Dune::XT::Grid::is_conforming_alugrid<Grid>::value ? (dim == 3 ? 1.1 : 1.06) : (dim == 3 ? 2.05e-1 : 1.45e-1);
  return tolerance;
}


#endif // DUNE_GDT_TEST_SPACES_RT_PDELAB_HH
