// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>
#include <vector>

#include <dune/common/typetraits.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/pdelab.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


// +---------------------------------------------------+
// | Definition of all Grids and GridParts of interest |
// +---------------------------------------------------+
typedef Dune::SGrid<1, 1> S1dGridType;
typedef Dune::grid::Part::Leaf::Const<S1dGridType> S1dGridPartType;
typedef Dune::SGrid<2, 2> S2dGridType;
typedef Dune::grid::Part::Leaf::Const<S2dGridType> S2dGridPartType;
typedef Dune::SGrid<3, 3> S3dGridType;
typedef Dune::grid::Part::Leaf::Const<S3dGridType> S3dGridPartType;

typedef Dune::YaspGrid<1> Yasp1dGridType;
typedef Dune::grid::Part::Leaf::Const<Yasp1dGridType> Yasp1dGridPartType;
typedef Dune::YaspGrid<2> Yasp2dGridType;
typedef Dune::grid::Part::Leaf::Const<Yasp2dGridType> Yasp2dGridPartType;
typedef Dune::YaspGrid<3> Yasp3dGridType;
typedef Dune::grid::Part::Leaf::Const<Yasp3dGridType> Yasp3dGridPartType;

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;
typedef Dune::grid::Part::Leaf::Const<AluConform2dGridType> AluConform2dGridPartType;
typedef Dune::ALUSimplexGrid<2, 2> AluSimplex2dGridType;
typedef Dune::grid::Part::Leaf::Const<AluSimplex2dGridType> AluSimplex2dGridPartType;
typedef Dune::ALUSimplexGrid<3, 3> AluSimplex3dGridType;
typedef Dune::grid::Part::Leaf::Const<AluSimplex3dGridType> AluSimplex3dGridPartType;
typedef Dune::ALUCubeGrid<3, 3> AluCube3dGridType;
typedef Dune::grid::Part::Leaf::Const<AluCube3dGridType> AluCube3dGridPartType;
#endif

typedef testing::
    Types<Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dGridPartType, 1, double, 3>
#if HAVE_ALUGRID
          ,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dGridPartType, 1, double, 3>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dGridPartType, 1, double, 2>,
          Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dGridPartType, 1, double, 3>
#endif
          ,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S1dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp1dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp3dGridPartType, 1, double, 1>
#if HAVE_ALUGRID
          ,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluConform2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluSimplex2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluSimplex3dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluCube3dGridPartType, 1, double, 1>
#endif
          ,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 1, double, 1>
#ifdef HAVE_ALUGRID
          ,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dGridPartType, 1, double, 1>,
          Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dGridPartType, 1, double, 1>
#endif
          ,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 1, double, 1>,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 1, double, 1>
#ifdef HAVE_ALUGRID
          ,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dGridPartType, 1, double, 1>,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dGridPartType, 1, double, 1>,
          Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dGridPartType, 1, double, 1>
#endif
          > SpaceTypes;

template <class SpaceType>
struct Spaces : public ::testing::Test
{
  typedef typename SpaceType::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;

  void crtp_check() const
  {
    // grid
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid_ptr  = grid_provider.grid();
    const auto grid_part = std::make_shared<const GridPartType>(*grid_ptr);
    // check the space
    const SpaceType space(grid_part);
    // check for static information
    typedef typename SpaceType::Traits Traits;
    typedef typename SpaceType::GridPartType S_GridPartType;
    typedef typename SpaceType::DomainFieldType S_DomainFieldType;
    static const unsigned int s_dimDomain = SpaceType::dimDomain;
    typedef typename SpaceType::RangeFieldType S_RangeFieldType;
    static const unsigned int s_dimRange     = SpaceType::dimRange;
    static const unsigned int s_dimRangeCols = SpaceType::dimRangeCols;
    static const int DUNE_UNUSED(s_polOrder) = SpaceType::polOrder;
    typedef typename SpaceType::BackendType S_BackendType;
    typedef typename SpaceType::MapperType MapperType;
    typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename SpaceType::EntityType EntityType;
    typedef typename SpaceType::PatternType PatternType;
    // check for functionality
    const auto entityIt      = grid_part->template begin<0>();
    const EntityType& entity = *entityIt;
    typedef typename Dune::GDT::SpaceInterface<Traits> SpaceInterfaceType;
    const SpaceInterfaceType& spaceAsInterface = static_cast<const SpaceInterfaceType&>(space);
    const std::shared_ptr<const S_GridPartType> DUNE_UNUSED(s_gridPart) = spaceAsInterface.gridPart();
    const S_BackendType& DUNE_UNUSED(s_backend) = spaceAsInterface.backend();
    const bool DUNE_UNUSED(s_continuous) = spaceAsInterface.continuous();
    const MapperType& mapper                  = spaceAsInterface.mapper();
    const BaseFunctionSetType baseFunctionSet = spaceAsInterface.baseFunctionSet(entity);
    const std::unique_ptr<const PatternType> pattern(spaceAsInterface.computePattern());

    // check the mapper for static information
    typedef typename MapperType::Traits M_Traits;
    typedef typename M_Traits::BackendType M_BackendType;
    // check the mapper for functionality
    typedef Dune::GDT::MapperInterface<M_Traits> MapperInterfaceType;
    const MapperInterfaceType& mapperAsInterface = static_cast<const MapperInterfaceType&>(mapper);
    const M_BackendType& DUNE_UNUSED(m_backend) = mapperAsInterface.backend();
    const size_t m_maxNumDofs = mapperAsInterface.maxNumDofs();
    const size_t DUNE_UNUSED(m_numDofs) = mapperAsInterface.numDofs(entity);
    Dune::DynamicVector<size_t> globalIndices(m_maxNumDofs, size_t(0));
    mapperAsInterface.globalIndices(entity, globalIndices);
    const size_t DUNE_UNUSED(globalIndex) = mapperAsInterface.mapToGlobal(entity, 0);

    // check the basefunctionset for static information
    typedef typename BaseFunctionSetType::Traits B_Traits;
    typedef typename B_Traits::BackendType B_BackendType;
    typedef typename B_Traits::EntityType B_EntityType;
    typedef typename BaseFunctionSetType::DomainType B_DomainType;
    typedef typename BaseFunctionSetType::RangeType B_RangeType;
    typedef typename BaseFunctionSetType::JacobianRangeType B_JacobianRangeType;
    // check the basefunctionset for functionality
    typedef Dune::GDT::
        BaseFunctionSetInterface<B_Traits, S_DomainFieldType, s_dimDomain, S_RangeFieldType, s_dimRange, s_dimRangeCols>
            BaseFunctionSetInterfaceType;
    const BaseFunctionSetInterfaceType& baseFunctionSetAsInterface =
        static_cast<const BaseFunctionSetInterfaceType&>(baseFunctionSet);
    const B_EntityType& DUNE_UNUSED(b_entity) = baseFunctionSetAsInterface.entity();
    const B_BackendType& DUNE_UNUSED(b_backend) = baseFunctionSetAsInterface.backend();
    const size_t b_size = baseFunctionSetAsInterface.size();
    const size_t DUNE_UNUSED(b_order) = baseFunctionSetAsInterface.order();
    const B_DomainType x = entity.geometry().center();
    std::vector<B_RangeType> values(b_size, B_RangeType(0));
    std::vector<B_JacobianRangeType> jacobians(b_size, B_JacobianRangeType(0));
    baseFunctionSetAsInterface.evaluate(x, values);
    baseFunctionSetAsInterface.jacobian(x, jacobians);
  } // ... crtp_check()
}; // struct Spaces


TYPED_TEST_CASE(Spaces, SpaceTypes);
TYPED_TEST(Spaces, crtp_check)
{
  this->crtp_check();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
