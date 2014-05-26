// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

// Then this one (otherwise we get alugrid problems)!
#include "spaces.hh"

#if !HAVE_DUNE_FEM && !HAVE_DUNE_PDELAB
#error "These tests requires at least one discretization module!"
#endif

#include <memory>
#include <vector>
#include <sstream>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/print.hh>

#include <dune/gdt/spaces/discontinuouslagrange/fem.hh>
#include <dune/gdt/spaces/discontinuouslagrange/pdelab.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


template <class SpaceType>
class Any_Space : public ::testing::Test, public SpaceTestBase<SpaceType>
{
};


#if HAVE_DUNE_FEM
#define DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                              \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 1, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 2, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 2, double, 1>

#if HAVE_ALUGRID
#define DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                      \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 2, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 2, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 2, double, 1>
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM

//#if HAVE_DUNE_PDELAB
//# define DISCONTINUOUS_LAGRANGE_SPACES_PDELAB \
//    Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< S1dLeafGridPartType, 1, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< S1dLeafGridPartType, 2, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< Yasp1dLeafGridPartType, 1, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< Yasp1dLeafGridPartType, 2, double, 1 >

//# if HAVE_ALUGRID
//#   define DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB \
//    Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluConform2dLeafGridPartType, 1, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluConform2dLeafGridPartType, 2, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluSimplex2dLeafGridPartType, 2, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluSimplex3dLeafGridPartType, 1, double, 1 > \
//  , Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased< AluSimplex3dLeafGridPartType, 2, double, 1 >
//# endif // HAVE_ALUGRID
//#endif //HAVE_DUNE_PDELAB


typedef testing::Types<DISCONTINUOUS_LAGRANGE_SPACES_FEM
//#if HAVE_DUNE_PDELAB && HAVE_DUNE_FEM
//                      ,
//#endif
//                      DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID && HAVE_DUNE_FEM
                       ,
#endif
                       DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
                       //#if HAVE_ALUGRID && HAVE_DUNE_PDELAB
                       //                      ,
                       //#endif
                       //                      DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
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
