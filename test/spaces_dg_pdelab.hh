// This file is View of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_DG_PDELAB_HH
#define DUNE_GDT_TEST_SPACES_DG_PDELAB_HH

#include <dune/gdt/playground/spaces/discontinuouslagrange/pdelab.hh>

#include "grids.hh"

#if HAVE_DUNE_PDELAB


#define SPACE_DG_PDELAB_SGRID(dd, rr, pp)                                                                              \
  Spaces::DiscontinuousLagrange::PdelabBased<S##dd##dLeafGridViewType, pp, double, rr>

#define SPACE_DG_PDELAB_YASPGRID(dd, rr, pp)                                                                           \
  Spaces::DiscontinuousLagrange::PdelabBased<Yasp##dd##dLeafGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB(pp)                                                                                           \
  SPACE_DG_PDELAB_SGRID(1, 1, pp)                                                                                      \
  , SPACE_DG_PDELAB_SGRID(2, 1, pp), SPACE_DG_PDELAB_SGRID(3, 1, pp), SPACE_DG_PDELAB_YASPGRID(1, 1, pp),              \
      SPACE_DG_PDELAB_YASPGRID(2, 1, pp), SPACE_DG_PDELAB_YASPGRID(3, 1, pp)


#if HAVE_ALUGRID


#define SPACE_DG_PDELAB_ALUCUBEGRID(dd, rr, pp)                                                                        \
  Spaces::DiscontinuousLagrange::PdelabBased<AluCube##dd##dLeafGridViewType, pp, double, rr>

#define SPACES_DG_PDELAB_ALUGRID(pp)                                                                                   \
  SPACE_DG_PDELAB_ALUCUBEGRID(2, 1, pp)                                                                                \
  , SPACE_DG_PDELAB_ALUCUBEGRID(3, 1, pp)


#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB

#endif // DUNE_GDT_TEST_SPACES_DG_PDELAB_HH
