// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_CG_PDELAB_HH
#define DUNE_GDT_TEST_SPACES_CG_PDELAB_HH

#include <dune/gdt/spaces/continuouslagrange/pdelab.hh>

#include "grids.hh"

#if HAVE_DUNE_PDELAB


#define SPACE_CG_PDELAB_SGRID(dd, rr, pp) \
  Spaces::ContinuousLagrange::PdelabBased< S ## dd ## dLeafGridViewType, pp, double, rr >

#define SPACE_CG_PDELAB_YASPGRID(dd, rr, pp) \
  Spaces::ContinuousLagrange::PdelabBased< Yasp ## dd ## dLeafGridViewType, pp, double, rr >

#define SPACES_CG_PDELAB(pp) \
    SPACE_CG_PDELAB_SGRID(1, pp, 1) \
  , SPACE_CG_PDELAB_SGRID(2, 1, pp) \
  , SPACE_CG_PDELAB_SGRID(3, 1, pp) \
  , SPACE_CG_PDELAB_YASPGRID(1, pp, 1) \
  , SPACE_CG_PDELAB_YASPGRID(2, 1, pp) \
  , SPACE_CG_PDELAB_YASPGRID(3, 1, pp)


# if HAVE_ALUGRID


#define SPACE_CG_PDELAB_ALUCONFORMGRID(dd, rr, pp) \
  Spaces::ContinuousLagrange::PdelabBased< AluConform ## dd ## dLeafGridViewType, pp, double, rr >

#define SPACE_CG_PDELAB_ALUCUBEGRID(dd, rr, pp) \
  Spaces::ContinuousLagrange::PdelabBased< AluCube ## dd ## dLeafGridViewType, pp, double, rr >

#define SPACES_CG_PDELAB_ALUGRID(pp) \
    SPACE_CG_PDELAB_ALUCONFORMGRID(2, 1, pp) \
  , SPACE_CG_PDELAB_ALUCONFORMGRID(3, 1, pp) \
  , SPACE_CG_PDELAB_ALUCUBEGRID(2, 1, pp) \
  , SPACE_CG_PDELAB_ALUCUBEGRID(3, 1, pp)


# endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB

#endif // DUNE_GDT_TEST_SPACES_CG_PDELAB_HH

