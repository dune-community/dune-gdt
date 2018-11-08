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

#ifndef DUNE_GDT_TEST_SPACES_CG_DEFAULT_HH
#define DUNE_GDT_TEST_SPACES_CG_DEFAULT_HH

// TODO: fix this test
#if 0
#include <dune/xt/grid/type_traits.hh>
#include <dune/gdt/test/grids.hh>
#include <dune/gdt/spaces/cg.hh>


#define SPACE_CG_YASPGRID(dd, rr, pp) Dune::GDT::ContinuousLagrangeSpace<Yasp##dd##dLeafGridViewType, pp, double>

#define SPACE_CG_YASPGRID_LEVEL(dd, rr, pp) Dune::GDT::ContinuousLagrangeSpace<Yasp##dd##dLevelGridViewType, pp, double>

#define SPACES_CG_LEVEL(pp)                                                                                            \
  SPACE_CG_YASPGRID_LEVEL(1, 1, pp)                                                                                    \
  , SPACE_CG_YASPGRID_LEVEL(2, 1, pp), SPACE_CG_YASPGRID_LEVEL(3, 1, pp)


#if HAVE_DUNE_ALUGRID


#define SPACE_CG_ALUCONFORMGRID(dd, rr, pp)                                                                            \
  Dune::GDT::ContinuousLagrangeSpace<AluConform##dd##dLeafGridViewType, pp, double>

#define SPACE_CG_ALUCUBEGRID(dd, rr, pp) Dune::GDT::ContinuousLagrangeSpace<AluCube##dd##dLeafGridViewType, pp, double>

#define SPACES_CG_ALUGRID(pp)                                                                                          \
  SPACE_CG_ALUCONFORMGRID(2, 1, pp)                                                                                    \
  , SPACE_CG_ALUCONFORMGRID(3, 1, pp), SPACE_CG_ALUCUBEGRID(2, 1, pp), SPACE_CG_ALUCUBEGRID(3, 1, pp)


#define SPACE_CG_ALUCUBEGRID_LEVEL(dd, rr, pp)                                                                         \
  Dune::GDT::ContinuousLagrangeSpace<AluCube##dd##dLevelGridViewType, pp, double>

#define SPACES_CG_ALUGRID_LEVEL(pp)                                                                                    \
  SPACE_CG_ALUCUBEGRID_LEVEL(2, 1, pp)                                                                                 \
  , SPACE_CG_ALUCUBEGRID_LEVEL(3, 1, pp)


#endif // HAVE_DUNE_ALUGRID

#endif // 0

#endif // DUNE_GDT_TEST_SPACES_CG_DEFAULT_HH
