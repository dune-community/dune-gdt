// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

// Then this one (otherwise we get alugrid problems)!
#include "spaces.hh"

#if !HAVE_DUNE_FEM && !HAVE_DUNE_FEM_LOCALFUNCTIONS && !HAVE_DUNE_PDELAB
#error "These tests requires at least one discretization module!"
#endif

#include <memory>
#include <vector>
#include <sstream>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/print.hh>

#include <dune/gdt/spaces/continuouslagrange/fem.hh>
#include <dune/gdt/spaces/continuouslagrange/pdelab.hh>
#include <dune/gdt/spaces/discontinuouslagrange/fem.hh>
#include <dune/gdt/spaces/discontinuouslagrange/pdelab.hh>
#include <dune/gdt/playground/spaces/raviartthomas/fem-localfunctions.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


template <class SpaceType>
class Any_Space : public ::testing::Test, public SpaceTestBase<SpaceType>
{
};


// typedef Dune::SGrid< 1, 1 > S1dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S1dGridType >::Type         S1dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S1dGridType, false >::Type  S1dGridPartType;
// typedef Dune::SGrid< 2, 2 > S2dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S2dGridType >::Type         S2dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S2dGridType, false >::Type  S2dGridPartType;
// typedef Dune::SGrid< 3, 3 > S3dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S3dGridType >::Type         S3dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< S3dGridType, false >::Type  S3dGridPartType;

// typedef Dune::YaspGrid< 1 > Yasp1dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp1dGridType >::Type         Yasp1dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp1dGridType, false >::Type  Yasp1dGridPartType;
// typedef Dune::YaspGrid< 2 > Yasp2dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp2dGridType >::Type         Yasp2dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp2dGridType, false >::Type  Yasp2dGridPartType;
// typedef Dune::YaspGrid< 3 > Yasp3dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp3dGridType >::Type         Yasp3dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< Yasp3dGridType, false >::Type  Yasp3dGridPartType;

//#if HAVE_ALUGRID
// typedef Dune::ALUConformGrid< 2, 2 > AluConform2dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType >::Type AluConform2dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, false >::Type
// AluConform2dGridPartType;
// typedef Dune::ALUSimplexGrid< 2, 2 > AluSimplex2dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex2dGridType >::Type AluSimplex2dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex2dGridType, false >::Type
// AluSimplex2dGridPartType;
// typedef Dune::ALUSimplexGrid< 3, 3 > AluSimplex3dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex3dGridType >::Type AluSimplex3dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex3dGridType, false >::Type
// AluSimplex3dGridPartType;
// typedef Dune::ALUCubeGrid< 3, 3 > AluCube3dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluCube3dGridType >::Type         AluCube3dGridViewType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluCube3dGridType, false >::Type  AluCube3dGridPartType;
//#endif

#if HAVE_ALUGRID
#define RTN0_RAVIART_THOMAS_SPACES_ALUGRID_FEM_LOCALFUNCTIONS                                                          \
  Dune::GDT::Spaces::RaviartThomas::FemLocalfunctionsBased<AluConform2dLeafGridPartType, 0, double, 2>,                \
      Dune::GDT::Spaces::RaviartThomas::FemLocalfunctionsBased<AluSimplex2dLeafGridPartType, 0, double, 2>
#endif // HAVE_ALUGRID

typedef testing::Types<

    > All_Spaces;


TYPED_TEST_CASE(Any_Space, All_Spaces);
TYPED_TEST(Any_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(Any_Space, All_Spaces);
TYPED_TEST(Any_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(Any_Space, All_Spaces);
TYPED_TEST(Any_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error:\n" << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
