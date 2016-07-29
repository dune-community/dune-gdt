// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_SPACES_DG_FEM_HH
#define DUNE_GDT_TEST_SPACES_DG_FEM_HH

#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>

#include <dune/gdt/test/grids.hh>

#if HAVE_DUNE_FEM


#define SPACE_DG_FEM_SGRID(dd, rr, pp) Dune::GDT::DuneFemDgSpaceWrapper<S##dd##dLeafGridPartType, pp, double, rr>

#define SPACE_DG_FEM_YASPGRID(dd, rr, pp) Dune::GDT::DuneFemDgSpaceWrapper<Yasp##dd##dLeafGridPartType, pp, double, rr>

#define SPACES_DG_FEM(pp)                                                                                              \
  SPACE_DG_FEM_SGRID(1, 1, pp)                                                                                         \
  , SPACE_DG_FEM_SGRID(2, 1, pp), SPACE_DG_FEM_SGRID(3, 1, pp), SPACE_DG_FEM_YASPGRID(1, 1, pp),                       \
      SPACE_DG_FEM_YASPGRID(2, 1, pp), SPACE_DG_FEM_YASPGRID(3, 1, pp)


#define SPACE_DG_FEM_SGRID_LEVEL(dd, rr, pp) Dune::GDT::DuneFemDgSpaceWrapper<S##dd##dLevelGridPartType, pp, double, rr>

#define SPACE_DG_FEM_YASPGRID_LEVEL(dd, rr, pp)                                                                        \
  Dune::GDT::DuneFemDgSpaceWrapper<Yasp##dd##dLevelGridPartType, pp, double, rr>

#define SPACES_DG_FEM_LEVEL(pp)                                                                                        \
  SPACE_DG_FEM_SGRID_LEVEL(1, 1, pp)                                                                                   \
  , SPACE_DG_FEM_SGRID_LEVEL(2, 1, pp), SPACE_DG_FEM_SGRID_LEVEL(3, 1, pp), SPACE_DG_FEM_YASPGRID_LEVEL(1, 1, pp),     \
      SPACE_DG_FEM_YASPGRID_LEVEL(2, 1, pp), SPACE_DG_FEM_YASPGRID_LEVEL(3, 1, pp)


#if HAVE_ALUGRID


#define SPACE_DG_FEM_ALUCONFORMGRID(dd, rr, pp)                                                                        \
  Dune::GDT::DuneFemDgSpaceWrapper<AluConform##dd##dLeafGridPartType, pp, double, rr>

#define SPACE_DG_FEM_ALUCUBEGRID(dd, rr, pp)                                                                           \
  Dune::GDT::DuneFemDgSpaceWrapper<AluCube##dd##dLeafGridPartType, pp, double, rr>

#define SPACES_DG_FEM_ALUGRID(pp)                                                                                      \
  SPACE_DG_FEM_ALUCONFORMGRID(2, 1, pp)                                                                                \
  , SPACE_DG_FEM_ALUCONFORMGRID(3, 1, pp), SPACE_DG_FEM_ALUCUBEGRID(2, 1, pp), SPACE_DG_FEM_ALUCUBEGRID(3, 1, pp)


#define SPACE_DG_FEM_ALUCONFORMGRID_LEVEL(dd, rr, pp)                                                                  \
  Dune::GDT::DuneFemDgSpaceWrapper<AluConform##dd##dLevelGridPartType, pp, double, rr>

#define SPACE_DG_FEM_ALUCUBEGRID_LEVEL(dd, rr, pp)                                                                     \
  Dune::GDT::DuneFemDgSpaceWrapper<AluCube##dd##dLevelGridPartType, pp, double, rr>

#define SPACES_DG_FEM_ALUGRID_LEVEL(pp)                                                                                \
  SPACE_DG_FEM_ALUCONFORMGRID_LEVEL(2, 1, pp)                                                                          \
  , SPACE_DG_FEM_ALUCONFORMGRID_LEVEL(3, 1, pp), SPACE_DG_FEM_ALUCUBEGRID_LEVEL(2, 1, pp),                             \
      SPACE_DG_FEM_ALUCUBEGRID_LEVEL(3, 1, pp)


#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM

#endif // DUNE_GDT_TEST_SPACES_DG_FEM_HH
