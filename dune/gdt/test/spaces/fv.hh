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

#ifndef DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
#define DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH

// TODO: reenable this test
#if 0
#  include <dune/gdt/spaces/l2/finite-volume.hh>

#  include <dune/gdt/test/grids.hh>


#  define SPACE_FV_YASPGRID(dd, rr) Dune::GDT::FvSpace<Yasp##dd##dLeafGridViewType, double, rr>

#  define SPACES_FV                                                                                                    \
    SPACE_FV_YASPGRID(1, 1)                                                                                            \
    , SPACE_FV_YASPGRID(1, 2), SPACE_FV_YASPGRID(1, 3), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(2, 2),              \
        SPACE_FV_YASPGRID(2, 3), SPACE_FV_YASPGRID(3, 1), SPACE_FV_YASPGRID(3, 2), SPACE_FV_YASPGRID(3, 3)

#  define SPACE_FV_YASPGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<Yasp##dd##dLevelGridViewType, double, rr>

#  if HAVE_DUNE_ALUGRID

#    define SPACE_FV_ALUCONFORMGRID(dd, rr) Dune::GDT::FvSpace<AluConform##dd##dLeafGridViewType, double, rr>

#    define SPACE_FV_ALUCUBEGRID(dd, rr) Dune::GDT::FvSpace<AluCube##dd##dLeafGridViewType, double, rr>

#    define SPACES_FV_ALUGRID                                                                                          \
      SPACE_FV_ALUCONFORMGRID(2, 1)                                                                                    \
      , SPACE_FV_ALUCONFORMGRID(2, 2), SPACE_FV_ALUCONFORMGRID(2, 3), SPACE_FV_ALUCONFORMGRID(3, 1),                   \
          SPACE_FV_ALUCONFORMGRID(3, 2), SPACE_FV_ALUCONFORMGRID(3, 3), SPACE_FV_ALUCUBEGRID(2, 1),                    \
          SPACE_FV_ALUCUBEGRID(2, 2), SPACE_FV_ALUCUBEGRID(2, 3), SPACE_FV_ALUCUBEGRID(3, 1),                          \
          SPACE_FV_ALUCUBEGRID(3, 2), SPACE_FV_ALUCUBEGRID(3, 3)


#    define SPACE_FV_ALUCONFORMGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<AluConform##dd##dLevelGridViewType, double, rr>

#    define SPACE_FV_ALUCUBEGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<AluCube##dd##dLevelGridViewType, double, rr>

#  endif // HAVE_DUNE_ALUGRID


#endif // 0
#endif // DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
