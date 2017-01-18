// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_SPACES_DG_PDELAB_HH
#define DUNE_GDT_TEST_SPACES_DG_PDELAB_HH

#include <dune/gdt/playground/spaces/dg/dune-pdelab-wrapper.hh>

#include <dune/gdt/test/grids.hh>

#if HAVE_DUNE_PDELAB


#define SPACE_DG_PDELAB_YASPGRID(dd, rr, pp)                                                                           \
  Dune::GDT::DunePdelabDgSpaceWrapper<Yasp##dd##dLeafGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB(pp)                                                                                           \
  SPACE_DG_PDELAB_YASPGRID(1, 1, pp), SPACE_DG_PDELAB_YASPGRID(2, 1, pp), SPACE_DG_PDELAB_YASPGRID(3, 1, pp)

#define SPACE_DG_PDELAB_YASPGRID_LEVEL(dd, rr, pp)                                                                     \
  Dune::GDT::DunePdelabDgSpaceWrapper<Yasp##dd##dLevelGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB_LEVEL(pp)                                                                                     \
  SPACE_DG_PDELAB_YASPGRID_LEVEL(1, 1, pp)                                                                             \
  , SPACE_DG_PDELAB_YASPGRID_LEVEL(2, 1, pp), SPACE_DG_PDELAB_YASPGRID_LEVEL(3, 1, pp)


#if HAVE_ALUGRID


#define SPACE_DG_PDELAB_ALUCUBEGRID(dd, rr, pp)                                                                        \
  Dune::GDT::DunePdelabDgSpaceWrapper<AluCube##dd##dLeafGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB_ALUGRID(pp)                                                                                   \
  SPACE_DG_PDELAB_ALUCUBEGRID(2, 1, pp)                                                                                \
  , SPACE_DG_PDELAB_ALUCUBEGRID(3, 1, pp)


#define SPACE_DG_PDELAB_ALUCUBEGRID_LEVEL(dd, rr, pp)                                                                  \
  Dune::GDT::DunePdelabDgSpaceWrapper<AluCube##dd##dLevelGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB_ALUGRID_LEVEL(pp)                                                                             \
  SPACE_DG_PDELAB_ALUCUBEGRID_LEVEL(2, 1, pp)                                                                          \
  , SPACE_DG_PDELAB_ALUCUBEGRID_LEVEL(3, 1, pp)


#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB

#endif // DUNE_GDT_TEST_SPACES_DG_PDELAB_HH
