// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_SPACES_DG_FEM_HH
#define DUNE_GDT_TEST_SPACES_DG_FEM_HH

#include <dune/gdt/playground/spaces/discontinuouslagrange/fem.hh>

#include "grids.hh"

#if HAVE_DUNE_FEM


#define SPACE_DG_FEM_SGRID(dd, rr, pp)                                                                                 \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S##dd##dLeafGridPartType, pp, double, rr>

#define SPACE_DG_FEM_YASPGRID(dd, rr, pp)                                                                              \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp##dd##dLeafGridPartType, pp, double, rr>

#define SPACES_DG_FEM(pp)                                                                                              \
  SPACE_DG_FEM_SGRID(1, 1, pp)                                                                                         \
  , SPACE_DG_FEM_SGRID(2, 1, pp), SPACE_DG_FEM_SGRID(3, 1, pp), SPACE_DG_FEM_YASPGRID(1, 1, pp),                       \
      SPACE_DG_FEM_YASPGRID(2, 1, pp), SPACE_DG_FEM_YASPGRID(3, 1, pp)


#if HAVE_ALUGRID


#define SPACE_DG_FEM_ALUSIMPLEXGRID(dd, rr, pp)                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex##dd##dLeafGridPartType, pp, double, rr>

#define SPACE_DG_FEM_ALUCUBEGRID(dd, rr, pp)                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube##dd##dLeafGridPartType, pp, double, rr>

#define SPACES_DG_FEM_ALUGRID(pp)                                                                                      \
  SPACE_DG_FEM_ALUSIMPLEXGRID(2, 1, pp)                                                                                \
  , SPACE_DG_FEM_ALUSIMPLEXGRID(3, 1, pp), SPACE_DG_FEM_ALUCUBEGRID(2, 1, pp), SPACE_DG_FEM_ALUCUBEGRID(3, 1, pp)


#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM

#endif // DUNE_GDT_TEST_SPACES_DG_FEM_HH
