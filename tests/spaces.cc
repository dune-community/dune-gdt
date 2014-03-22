// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>
#include <vector>
#include <sstream>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/print.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/pdelab.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


template< class SpaceType >
struct Any_Space
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;

  void fulfills_interface() const
  {
    // grid
    const GridProviderType grid_provider(0.0, 1.0, 2u);
    const auto grid_ptr = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid_ptr);
    // check the space
    const SpaceType space(grid_part);
    // check for static information
    typedef typename SpaceType::Traits              Traits;
    typedef typename SpaceType::GridPartType        S_GridPartType;
    typedef typename SpaceType::DomainFieldType     S_DomainFieldType;
    static const unsigned int                       s_dimDomain = SpaceType::dimDomain;
    typedef typename SpaceType::RangeFieldType      S_RangeFieldType;
    static const unsigned int                       s_dimRange = SpaceType::dimRange;
    static const unsigned int                       s_dimRangeCols = SpaceType::dimRangeCols;
    static const int DUNE_UNUSED(                   s_polOrder) = SpaceType::polOrder;
    typedef typename SpaceType::BackendType         S_BackendType;
    typedef typename SpaceType::MapperType          MapperType;
    typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename SpaceType::EntityType          EntityType;
    typedef typename SpaceType::PatternType         PatternType;
    // check for functionality
    const auto entityIt = grid_part->template begin< 0 >();
    const EntityType& entity = *entityIt;
    typedef typename Dune::GDT::SpaceInterface< Traits > SpaceInterfaceType;
    const SpaceInterfaceType& spaceAsInterface = static_cast< const SpaceInterfaceType& >(space);
    const std::shared_ptr< const S_GridPartType > DUNE_UNUSED(s_gridPart) = spaceAsInterface.gridPart();
    const S_BackendType& DUNE_UNUSED(       s_backend)      = spaceAsInterface.backend();
    const bool DUNE_UNUSED(                 s_continuous)   = spaceAsInterface.continuous();
    const MapperType&                       mapper          = spaceAsInterface.mapper();
    const BaseFunctionSetType               baseFunctionSet = spaceAsInterface.baseFunctionSet(entity);
    const std::unique_ptr< const PatternType > pattern(spaceAsInterface.computePattern());

    // check the mapper for static information
    typedef typename MapperType::Traits M_Traits;
    typedef typename M_Traits::BackendType M_BackendType;
    // check the mapper for functionality
    typedef Dune::GDT::MapperInterface< M_Traits > MapperInterfaceType;
    const MapperInterfaceType& mapperAsInterface = static_cast< const MapperInterfaceType& >(mapper);
    const M_BackendType& DUNE_UNUSED( m_backend)    = mapperAsInterface.backend();
    const size_t                      m_maxNumDofs  = mapperAsInterface.maxNumDofs();
    const size_t DUNE_UNUSED(         m_numDofs)    = mapperAsInterface.numDofs(entity);
    Dune::DynamicVector< size_t > globalIndices(m_maxNumDofs, size_t(0));
    mapperAsInterface.globalIndices(entity, globalIndices);
    const size_t DUNE_UNUSED(         globalIndex)  = mapperAsInterface.mapToGlobal(entity, 0);

    // check the basefunctionset for static information
    typedef typename BaseFunctionSetType::Traits  B_Traits;
    typedef typename B_Traits::BackendType        B_BackendType;
    typedef typename B_Traits::EntityType         B_EntityType;
    typedef typename BaseFunctionSetType::DomainType         B_DomainType;
    typedef typename BaseFunctionSetType::RangeType          B_RangeType;
    typedef typename BaseFunctionSetType::JacobianRangeType  B_JacobianRangeType;
    // check the basefunctionset for functionality
    typedef Dune::GDT::BaseFunctionSetInterface< B_Traits,
                                          S_DomainFieldType, s_dimDomain,
                                          S_RangeFieldType, s_dimRange, s_dimRangeCols > BaseFunctionSetInterfaceType;
    const BaseFunctionSetInterfaceType& baseFunctionSetAsInterface
        = static_cast< const BaseFunctionSetInterfaceType& >(baseFunctionSet);
    const B_EntityType& DUNE_UNUSED(  b_entity)  = baseFunctionSetAsInterface.entity();
    const B_BackendType& DUNE_UNUSED( b_backend) = baseFunctionSetAsInterface.backend();
    const size_t                      b_size     = baseFunctionSetAsInterface.size();
    const size_t       DUNE_UNUSED(   b_order)   = baseFunctionSetAsInterface.order();
    const B_DomainType x = entity.geometry().center();
    std::vector< B_RangeType > values(b_size, B_RangeType(0));
    std::vector< B_JacobianRangeType > jacobians(b_size, B_JacobianRangeType(0));
    baseFunctionSetAsInterface.evaluate(x, values);
    baseFunctionSetAsInterface.jacobian(x, jacobians);
  } // ... fulfills_interface()
}; // struct Any_Space


template< class SpaceType >
struct P1Q1_Continuous_Lagrange
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  static const unsigned int dimDomain = GridType::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  static std::vector< DomainFieldType > convert_vector(const DomainType& source)
  {
    std::vector< DomainFieldType > ret(dimDomain, DomainFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret[ii] = source[ii];
    return ret;
  }

  void maps_correctly()
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    grid_provider.visualize("grid");
    const auto grid_ptr = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid_ptr);
    const SpaceType space(grid_part);
    // walk the grid to create a map of all vertices
    std::map< std::vector< DomainFieldType >, std::set< size_t > > vertex_to_indices_map;
    const auto entity_end_it = grid_part->template end< 0 >();
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      for (int cc = 0; cc < entity.template count< dimDomain >(); ++cc) {
        const auto vertex_ptr = entity.template subEntity< dimDomain >(cc);
        const DomainType vertex = vertex_ptr->geometry().center();
        vertex_to_indices_map[convert_vector(vertex)] = std::set< size_t >();
      }
    }
    // walk the grid again to find all DoF ids
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      const int num_vertices = entity.template count< dimDomain >();
      const auto basis = space.baseFunctionSet(entity);
      if (basis.size() != num_vertices)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, "basis.size() = " << basis.size());
      for (int cc = 0; cc < num_vertices; ++cc) {
        const auto vertex_ptr = entity.template subEntity< dimDomain >(cc);
        const DomainType vertex = vertex_ptr->geometry().center();
        // find the local basis function which corresponds to this vertex
        const auto basis_values = basis.evaluate(entity.geometry().local(vertex));
        if (basis_values.size() != num_vertices)
          DUNE_THROW_COLORFULLY(Exceptions::internal_error, "basis_values.size() = " << basis_values.size());
        size_t matches = 0;
        size_t misses = 0;
        size_t local_DoF_index = 0;
        for (size_t ii = 0; ii < num_vertices; ++ii) {
          if (Common::FloatCmp::eq(basis_values[ii][0], typename SpaceType::RangeFieldType(1))) {
            local_DoF_index = ii;
            ++matches;
          } else
            ++misses;
        }
        if (matches != 1 || misses != (num_vertices - 1)) {
          std::stringstream ss;
          ss << "matches = " << matches << ", misses = " << misses << ", num_vertices = "
             << num_vertices << ", entity " << grid_part->indexSet().index(entity)
             << ", vertex " << cc << ": [ " << vertex << "], ";
          Common::print(basis_values, "basis_values", ss);
          DUNE_THROW_COLORFULLY(Exceptions::internal_error, ss.str());
        }
        // now we know that the local DoF index of this vertex is ii
        const size_t global_DoF_index = space.mapper().mapToGlobal(entity, local_DoF_index);
        vertex_to_indices_map[convert_vector(vertex)].insert(global_DoF_index);
      }
    }
    // check that all vertices have indeed one and only one global DoF id and that the numbering is consecutive
    std::set< size_t > global_DoF_indices;
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids = entry.second;
      if (vertex_ids.size() != 1)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, vertex_ids.size());
      global_DoF_indices.insert(*(vertex_ids.begin()));
    }
    if (vertex_to_indices_map.size() != global_DoF_indices.size())
      DUNE_THROW_COLORFULLY(Exceptions::internal_error,
                            "vertex_to_indices_map.size() = " << vertex_to_indices_map.size()
                            << ", global_DoF_indices.size() = " << global_DoF_indices.size());
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      if (global_DoF_id != count)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, "count = " << count << ", global_DoF_id = " << global_DoF_id);
      ++count;
    }
  } // ... maps_correctly()
}; // struct P1Q1_Continuous_Lagrange


typedef Dune::SGrid< 1, 1 >                           S1dGridType;
typedef Dune::grid::Part::Leaf::Const< S1dGridType >  S1dGridPartType;
typedef Dune::SGrid< 2, 2 >                           S2dGridType;
typedef Dune::grid::Part::Leaf::Const< S2dGridType >  S2dGridPartType;
typedef Dune::SGrid< 3, 3 >                           S3dGridType;
typedef Dune::grid::Part::Leaf::Const< S3dGridType >  S3dGridPartType;

typedef Dune::YaspGrid< 1 >                             Yasp1dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp1dGridType > Yasp1dGridPartType;
typedef Dune::YaspGrid< 2 >                             Yasp2dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp2dGridType > Yasp2dGridPartType;
typedef Dune::YaspGrid< 3 >                             Yasp3dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp3dGridType > Yasp3dGridPartType;

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid< 2, 2 >                          AluConform2dGridType;
typedef Dune::grid::Part::Leaf::Const< AluConform2dGridType > AluConform2dGridPartType;
typedef Dune::ALUSimplexGrid< 2, 2 >                          AluSimplex2dGridType;
typedef Dune::grid::Part::Leaf::Const< AluSimplex2dGridType > AluSimplex2dGridPartType;
typedef Dune::ALUSimplexGrid< 3, 3 >                          AluSimplex3dGridType;
typedef Dune::grid::Part::Leaf::Const< AluSimplex3dGridType > AluSimplex3dGridPartType;
typedef Dune::ALUCubeGrid< 3, 3 >                             AluCube3dGridType;
typedef Dune::grid::Part::Leaf::Const< AluCube3dGridType >    AluCube3dGridPartType;
#endif

#define P1_CONTINUOUS_LAGRANGE_SPACES \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< Yasp1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< Yasp1dGridPartType, 1, double, 1 >

#if HAVE_ALUGRID
# define P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< AluConform2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< AluSimplex2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< AluSimplex3dGridPartType, 1, double, 1 >
#endif // HAVE_ALUGRID

#define P2_CONTINUOUS_LAGRANGE_SPACES \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 2 >

#if HAVE_ALUGRID
# define P2_CONTINUOUS_LAGRANGE_SPACES_ALUGRID \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex2dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex3dGridPartType, 1, double, 2 >
#endif // HAVE_ALUGRID

#define Q1_CONTINUOUS_LAGRANGE_SPACES \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< Yasp1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< S1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< S2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< S3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< Yasp1dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< Yasp2dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< Yasp3dGridPartType, 1, double, 1 >

#if HAVE_ALUGRID
# define Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluCube3dGridPartType, 1, double, 1 > \
  , Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper< AluCube3dGridPartType, 1, double, 1 >
#endif // HAVE_ALUGRID

#define Q2_CONTINUOUS_LAGRANGE_SPACES \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 2 > \
  , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 2 >

#if HAVE_ALUGRID
# define Q2_CONTINUOUS_LAGRANGE_SPACES_ALUGRID \
    Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluCube3dGridPartType, 1, double, 2 >
#endif // HAVE_ALUGRID

typedef testing::Types< P1_CONTINUOUS_LAGRANGE_SPACES
#if HAVE_ALUGRID
                      , P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      , Q1_CONTINUOUS_LAGRANGE_SPACES
#if P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
                      , Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      > P1Q1_Continuous_Lagrange_Spaces;

typedef testing::Types< P1_CONTINUOUS_LAGRANGE_SPACES
#if HAVE_ALUGRID
                      , P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      , P2_CONTINUOUS_LAGRANGE_SPACES
#if HAVE_ALUGRID
                      , P2_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      , Q1_CONTINUOUS_LAGRANGE_SPACES
#if HAVE_ALUGRID
                      , Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      , Q2_CONTINUOUS_LAGRANGE_SPACES
#if HAVE_ALUGRID
                      , Q2_CONTINUOUS_LAGRANGE_SPACES_ALUGRID
#endif
                      > All_Spaces;

//                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< S1dGridPartType, 1, double, 1 >
//                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< Yasp1dGridPartType, 1, double, 1 >
//#if HAVE_ALUGRID
//                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
//                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 1, double, 1 >
//                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 1, double, 1 >
//#endif


TYPED_TEST_CASE(Any_Space, All_Spaces);
TYPED_TEST(Any_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, maps_correctly)
{
  this->maps_correctly();
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
} // ... main(...)
