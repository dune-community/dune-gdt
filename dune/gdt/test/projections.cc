// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include "projections.hh"

using namespace Dune::GDT::Test;

#if 0


 typedef testing::Types<SPACES_CG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;
 typedef testing::Types<SPACES_CG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_CG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;
 typedef testing::Types<SPACES_DG_FEM(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_DG_FEM_ALUGRID(1)
#endif
                       >
    SpaceTypes;
 typedef testing::Types<SPACES_DG_PDELAB(1)
#if HAVE_ALUGRID
                           ,
                       SPACES_DG_PDELAB_ALUGRID(1)
#endif
                       >
    SpaceTypes;
 typedef testing::Types<SPACE_FV_YASPGRID(1, 1), SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)
#if HAVE_ALUGRID
                                                                             ,
                       SPACE_FV_ALUCONFORMGRID(2, 1), SPACE_FV_ALUCONFORMGRID(3, 1), SPACE_FV_ALUCUBEGRID(2, 1),
                       SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_ALUGRID
                       >
    SpaceTypes;
#endif // 0

TEST_F(ProjectionTest, produces_correct_results)
{
  this->produces_correct_results();
}
