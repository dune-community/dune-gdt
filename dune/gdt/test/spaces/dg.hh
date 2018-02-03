// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_TEST_SPACES_DG_DEFAULT_HH
#define DUNE_GDT_TEST_SPACES_DG_DEFAULT_HH

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/test/grids.hh>


#define SPACE_DG_YASPGRID(dd, rr, pp) Dune::GDT::DiscontinuousLagrangeSpace<Yasp##dd##dLeafGridViewType, pp, double>

#define SPACE_DG_YASPGRID_LEVEL(dd, rr, pp)                                                                            \
  Dune::GDT::DiscontinuousLagrangeSpace<Yasp##dd##dLevelGridViewType, pp, double>

#define SPACES_DG_LEVEL(pp)                                                                                            \
  SPACE_DG_YASPGRID_LEVEL(1, 1, pp), SPACE_DG_YASPGRID_LEVEL(2, 1, pp), SPACE_DG_YASPGRID_LEVEL(3, 1, pp)


#if HAVE_DUNE_ALUGRID


#define SPACE_DG_ALUCONFORMGRID(dd, rr, pp)                                                                            \
  Dune::GDT::DiscontinuousLagrangeSpace<AluConform##dd##dLeafGridViewType, pp, double>

#define SPACE_DG_ALUCUBEGRID(dd, rr, pp)                                                                               \
  Dune::GDT::DiscontinuousLagrangeSpace<AluCube##dd##dLeafGridViewType, pp, double>

#define SPACES_DG_ALUGRID(pp)                                                                                          \
  SPACE_DG_ALUCONFORMGRID(2, 1, pp)                                                                                    \
  , SPACE_DG_ALUCONFORMGRID(3, 1, pp), SPACE_DG_ALUCUBEGRID(2, 1, pp), SPACE_DG_ALUCUBEGRID(3, 1, pp)


#define SPACE_DG_ALUCONFORMGRID_LEVEL(dd, rr, pp)                                                                      \
  Dune::GDT::DiscontinuousLagrangeSpace<AluConform##dd##dLevelGridViewType, pp, double>

#define SPACE_DG_ALUCUBEGRID_LEVEL(dd, rr, pp)                                                                         \
  Dune::GDT::DiscontinuousLagrangeSpace<AluCube##dd##dLevelGridViewType, pp, double>

#define SPACES_DG_ALUGRID_LEVEL(pp)                                                                                    \
  SPACE_DG_ALUCONFORMGRID_LEVEL(2, 1, pp)                                                                              \
  , SPACE_DG_ALUCONFORMGRID_LEVEL(3, 1, pp), SPACE_DG_ALUCUBEGRID_LEVEL(2, 1, pp), SPACE_DG_ALUCUBEGRID_LEVEL(3, 1, pp)


#endif // HAVE_DUNE_ALUGRID

#endif // DUNE_GDT_TEST_SPACES_DG_DEFAULT_HH
