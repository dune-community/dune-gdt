// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
#define DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH

#include <dune/gdt/spaces/fv/default.hh>

#include <dune/gdt/test/grids.hh>


#define SPACE_FV_SGRID(dd, rr) Dune::GDT::FvSpace<S##dd##dLeafGridViewType, double, rr>

#define SPACE_FV_YASPGRID(dd, rr) Dune::GDT::FvSpace<Yasp##dd##dLeafGridViewType, double, rr>

#define SPACES_FV                                                                                                      \
  SPACE_FV_SGRID(1, 1)                                                                                                 \
  , SPACE_FV_SGRID(1, 2), SPACE_FV_SGRID(1, 3), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(2, 2), SPACE_FV_SGRID(2, 3),      \
      SPACE_FV_SGRID(3, 1), SPACE_FV_SGRID(3, 2), SPACE_FV_SGRID(3, 2), SPACE_FV_YASPGRID(1, 1),                       \
      SPACE_FV_YASPGRID(1, 2), SPACE_FV_YASPGRID(1, 3), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(2, 2),              \
      SPACE_FV_YASPGRID(2, 3), SPACE_FV_YASPGRID(3, 1), SPACE_FV_YASPGRID(3, 2), SPACE_FV_YASPGRID(3, 3)


#define SPACE_FV_SGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<S##dd##dLevelGridViewType, double, rr>

#define SPACE_FV_YASPGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<Yasp##dd##dLevelGridViewType, double, rr>


#if HAVE_DUNE_ALUGRID

#define SPACE_FV_ALUCONFORMGRID(dd, rr) Dune::GDT::FvSpace<AluConform##dd##dLeafGridViewType, double, rr>

#define SPACE_FV_ALUCUBEGRID(dd, rr) Dune::GDT::FvSpace<AluCube##dd##dLeafGridViewType, double, rr>

#define SPACES_FV_ALUGRID                                                                                              \
  SPACE_FV_ALUCONFORMGRID(2, 1)                                                                                        \
  , SPACE_FV_ALUCONFORMGRID(2, 2), SPACE_FV_ALUCONFORMGRID(2, 3), SPACE_FV_ALUCUBEGRID(2, 1),                        \
      SPACE_FV_ALUCUBEGRID(2, 2), SPACE_FV_ALUCUBEGRID(2, 3), SPACE_FV_ALUCUBEGRID(3, 1), SPACE_FV_ALUCUBEGRID(3, 2),  \
      SPACE_FV_ALUCUBEGRID(3, 3)


#define SPACE_FV_ALUCONFORMGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<AluConform##dd##dLevelGridViewType, double, rr>

#define SPACE_FV_ALUCUBEGRID_LEVEL(dd, rr) Dune::GDT::FvSpace<AluCube##dd##dLevelGridViewType, double, rr>

#endif // HAVE_DUNE_ALUGRID


#endif // DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
