// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
#define DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH

#include "grids.hh"

#include <dune/gdt/playground/spaces/finitevolume/default.hh>


#define SPACE_FV_SGRID(dd, rr) Spaces::FiniteVolume::Default<S##dd##dLeafGridViewType, double, rr>

#define SPACE_FV_YASPGRID(dd, rr) Spaces::FiniteVolume::Default<Yasp##dd##dLeafGridViewType, double, rr>

#define SPACES_FV                                                                                                      \
  SPACE_FV_SGRID(1, 1)                                                                                                 \
  , SPACE_FV_SGRID(1, 2), SPACE_FV_SGRID(1, 3), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(2, 2), SPACE_FV_SGRID(2, 3),      \
      SPACE_FV_SGRID(3, 1), SPACE_FV_SGRID(3, 2), SPACE_FV_SGRID(3, 2), SPACE_FV_YASPGRID(1, 1),                       \
      SPACE_FV_YASPGRID(1, 2), SPACE_FV_YASPGRID(1, 3), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(2, 2),              \
      SPACE_FV_YASPGRID(2, 3), SPACE_FV_YASPGRID(3, 1), SPACE_FV_YASPGRID(3, 2), SPACE_FV_YASPGRID(3, 3)


#if HAVE_ALUGRID

#define SPACE_FV_ALUSIMPLEXGRID(dd, rr) Spaces::FiniteVolume::Default<AluSimplex##dd##dLeafGridViewType, double, rr>

#define SPACE_FV_ALUCUBEGRID(dd, rr) Spaces::FiniteVolume::Default<AluCube##dd##dLeafGridViewType, double, rr>

#define SPACES_FV_ALUGRID                                                                                              \
  SPACE_FV_ALUSIMPLEXGRID(2, 1)                                                                                        \
  , SPACE_FV_ALUSIMPLEXGRID(2, 2), SPACE_FV_ALUSIMPLEXGRID(2, 3), SPACE_FV_ALUSIMPLEXGRID(3, 1),                       \
      SPACE_FV_ALUSIMPLEXGRID(3, 2), SPACE_FV_ALUSIMPLEXGRID(3, 3), SPACE_FV_ALUCUBEGRID(2, 1),                        \
      SPACE_FV_ALUCUBEGRID(2, 2), SPACE_FV_ALUCUBEGRID(2, 3), SPACE_FV_ALUCUBEGRID(3, 1), SPACE_FV_ALUCUBEGRID(3, 2),  \
      SPACE_FV_ALUCUBEGRID(3, 3)

#endif // HAVE_ALUGRID


#endif // DUNE_GDT_TEST_SPACES_FV_DEFAULT_HH
