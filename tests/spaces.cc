// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#if !HAVE_DUNE_FEM && !HAVE_DUNE_FEM_LOCALFUNCTIONS && !HAVE_DUNE_PDELAB
#error "These tests requires at least one discretization module!"
#endif

#include <memory>
#include <vector>
#include <sstream>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/print.hh>

#include <dune/gdt/space/tools.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/pdelab.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


template <class SpaceType>
struct Any_Space : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;

  void fulfills_interface() const
  {
    // grid
    GridProviderType grid_provider(0.0, 1.0, 2u);
    auto grid_ptr             = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_ptr);
    // check the space
    const SpaceType space(grid_part_view);
    // check for static information
    typedef typename SpaceType::Traits Traits;
    typedef typename SpaceType::GridViewType S_GridViewType;
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
    const auto entityIt      = grid_part_view->template begin<0>();
    const EntityType& entity = *entityIt;
    typedef typename Dune::GDT::SpaceInterface<Traits> SpaceInterfaceType;
    const SpaceInterfaceType& spaceAsInterface = static_cast<const SpaceInterfaceType&>(space);
    const std::shared_ptr<const S_GridViewType> DUNE_UNUSED(s_gridView) = spaceAsInterface.grid_view();
    const S_BackendType& DUNE_UNUSED(s_backend) = spaceAsInterface.backend();
    const MapperType& mapper            = spaceAsInterface.mapper();
    BaseFunctionSetType baseFunctionSet = spaceAsInterface.base_function_set(entity);
    PatternType pattern                 = spaceAsInterface.compute_pattern();

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
  } // ... fulfills_interface()
}; // struct Any_Space


template <class SpaceType>
struct P1Q1_Continuous_Lagrange : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::GridProviderCube<GridType> GridProviderType;
  static const unsigned int dimDomain = GridType::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  static std::vector<DomainFieldType> convert_vector(const DomainType& source)
  {
    std::vector<DomainFieldType> ret(dimDomain, DomainFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret[ii] = source[ii];
    return ret;
  }

  template <class T, int d, class R, int r, int rC>
  void matches_signature(const Dune::GDT::ContinuousLagrangeSpaceBase<T, d, R, r, rC>& /*space*/)
  {
    static_assert(std::is_same<typename SpaceType::Traits, T>::value, "");
    static_assert(std::is_same<typename SpaceType::RangeFieldType, R>::value, "");
    static_assert(d == SpaceType::dimDomain, "");
    static_assert(r == SpaceType::dimRange, "");
    static_assert(rC == SpaceType::dimRangeCols, "");
  }

  void fulfills_interface()
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid_ptr             = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_ptr);
    const SpaceType space(grid_part_view);
    matches_signature(space);
    const auto entity_ptr                   = grid_part_view->template begin<0>();
    const auto& entity                      = *entity_ptr;
    const auto basis                        = space.base_function_set(entity);
    std::vector<DomainType> lagrange_points = space.lagrange_points(entity);
    if (lagrange_points.size() != basis.size())
      DUNE_THROW_COLORFULLY(
          Exceptions::internal_error,
          "lagrange_points.size() = " << lagrange_points.size() << ", basis.size() = " << basis.size());
    typedef typename SpaceType::IntersectionType IntersectionType;
    typedef typename SpaceType::RangeFieldType RangeFieldType;
    GridboundaryAllDirichlet<IntersectionType> boundary_info;
    std::set<size_t> local_dirichlet_DoFs = space.local_dirichlet_DoFs(entity, boundary_info);
    Constraints::Dirichlet<IntersectionType, RangeFieldType, true> dirichlet_constraints_a(boundary_info, 0, 0);
    Constraints::Dirichlet<IntersectionType, RangeFieldType, false> dirichlet_constraints_b(boundary_info, 0, 0);
    space.local_constraints(entity, dirichlet_constraints_a);
    space.local_constraints(entity, dirichlet_constraints_b);
  }

  void maps_correctly()
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto grid_ptr             = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_ptr);
    const SpaceType space(grid_part_view);
    // walk the grid to create a map of all vertices
    std::map<std::vector<DomainFieldType>, std::set<size_t>> vertex_to_indices_map;
    const auto entity_end_it = grid_part_view->template end<0>();
    for (auto entity_it = grid_part_view->template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      for (int cc = 0; cc < entity.template count<dimDomain>(); ++cc) {
        const auto vertex_ptr   = entity.template subEntity<dimDomain>(cc);
        const DomainType vertex = vertex_ptr->geometry().center();
        vertex_to_indices_map[convert_vector(vertex)] = std::set<size_t>();
      }
    }
    // walk the grid again to find all DoF ids
    for (auto entity_it = grid_part_view->template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity     = *entity_it;
      const int num_vertices = entity.template count<dimDomain>();
      const auto basis = space.base_function_set(entity);
      if (basis.size() != size_t(num_vertices))
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, "basis.size() = " << basis.size());
      for (int cc = 0; cc < num_vertices; ++cc) {
        const auto vertex_ptr   = entity.template subEntity<dimDomain>(cc);
        const DomainType vertex = vertex_ptr->geometry().center();
        // find the local basis function which corresponds to this vertex
        const auto basis_values = basis.evaluate(entity.geometry().local(vertex));
        if (basis_values.size() != size_t(num_vertices))
          DUNE_THROW_COLORFULLY(Exceptions::internal_error, "basis_values.size() = " << basis_values.size());
        size_t ones            = 0;
        size_t zeros           = 0;
        size_t failures        = 0;
        size_t local_DoF_index = 0;
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          if (Common::FloatCmp::eq(basis_values[ii][0], typename SpaceType::RangeFieldType(1))) {
            local_DoF_index = ii;
            ++ones;
          } else if (Common::FloatCmp::eq(basis_values[ii][0], typename SpaceType::RangeFieldType(0)))
            ++zeros;
          else
            ++failures;
        }
        if (ones != 1 || zeros != (basis.size() - 1) || failures > 0) {
          std::stringstream ss;
          ss << "ones = " << ones << ", zeros = " << zeros << ", failures = " << failures
             << ", num_vertices = " << num_vertices << ", entity " << grid_part_view->indexSet().index(entity)
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
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids = entry.second;
      if (vertex_ids.size() != 1)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, vertex_ids.size());
      global_DoF_indices.insert(*(vertex_ids.begin()));
    }
    if (vertex_to_indices_map.size() != global_DoF_indices.size())
      DUNE_THROW_COLORFULLY(Exceptions::internal_error,
                            "vertex_to_indices_map.size() = " << vertex_to_indices_map.size()
                                                              << ", global_DoF_indices.size() = "
                                                              << global_DoF_indices.size());
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      if (global_DoF_id != count)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, "count = " << count << ", global_DoF_id = " << global_DoF_id);
      ++count;
    }
  } // ... maps_correctly()
}; // struct P1Q1_Continuous_Lagrange


typedef Dune::SGrid<1, 1> S1dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S1dGridType>::Type S1dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S1dGridType, false>::Type S1dGridPartType;
typedef Dune::SGrid<2, 2> S2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S2dGridType>::Type S2dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S2dGridType, false>::Type S2dGridPartType;
typedef Dune::SGrid<3, 3> S3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S3dGridType>::Type S3dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S3dGridType, false>::Type S3dGridPartType;

typedef Dune::YaspGrid<1> Yasp1dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp1dGridType>::Type Yasp1dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp1dGridType, false>::Type Yasp1dGridPartType;
typedef Dune::YaspGrid<2> Yasp2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp2dGridType>::Type Yasp2dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp2dGridType, false>::Type Yasp2dGridPartType;
typedef Dune::YaspGrid<3> Yasp3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp3dGridType>::Type Yasp3dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp3dGridType, false>::Type Yasp3dGridPartType;

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType>::Type AluConform2dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dGridPartType;
typedef Dune::ALUSimplexGrid<2, 2> AluSimplex2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType>::Type AluSimplex2dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dGridPartType;
typedef Dune::ALUSimplexGrid<3, 3> AluSimplex3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType>::Type AluSimplex3dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dGridPartType;
typedef Dune::ALUCubeGrid<3, 3> AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType>::Type AluCube3dGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, false>::Type AluCube3dGridPartType;
#endif

#define P1_CONTINUOUS_LAGRANGE_SPACES_FEM                                                                              \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dGridPartType, 1, double, 1>,                                       \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dGridPartType, 1, double, 1>

#define P1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS                                                               \
  Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 1, double, 1>,                         \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 1, double, 1>

#define P1_CONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                           \
  Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S1dGridViewType, 1, double, 1>,                                    \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp1dGridViewType, 1, double, 1>

#if HAVE_ALUGRID
#define P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                      \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluConform2dGridPartType, 1, double, 1>,                              \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex2dGridPartType, 1, double, 1>,                          \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluSimplex3dGridPartType, 1, double, 1>

#define P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM_LOCALFUNCTIONS                                                       \
  Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dGridPartType, 1, double, 1>,                \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dGridPartType, 1, double, 1>,            \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex3dGridPartType, 1, double, 1>

#define P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                   \
  Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluConform2dGridViewType, 1, double, 1>,                           \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluSimplex2dGridViewType, 1, double, 1>,                       \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluSimplex3dGridViewType, 1, double, 1>
#endif // HAVE_ALUGRID

#define Q1_CONTINUOUS_LAGRANGE_SPACES_FEM                                                                              \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S1dGridPartType, 1, double, 1>,                                       \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S2dGridPartType, 1, double, 1>,                                   \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<S3dGridPartType, 1, double, 1>,                                   \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp1dGridPartType, 1, double, 1>,                                \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp2dGridPartType, 1, double, 1>,                                \
      Dune::GDT::ContinuousLagrangeSpace::FemWrapper<Yasp3dGridPartType, 1, double, 1>

#define Q1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS                                                               \
  Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 1, double, 1>,                         \
      Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 1, double, 1>

#define Q1_CONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                           \
  Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S1dGridViewType, 1, double, 1>,                                    \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S2dGridViewType, 1, double, 1>,                                \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<S3dGridViewType, 1, double, 1>,                                \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp1dGridViewType, 1, double, 1>,                             \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp2dGridViewType, 1, double, 1>,                             \
      Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<Yasp3dGridViewType, 1, double, 1>

#if HAVE_ALUGRID
#define Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                      \
  Dune::GDT::ContinuousLagrangeSpace::FemWrapper<AluCube3dGridPartType, 1, double, 1>

#define Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                   \
  Dune::GDT::ContinuousLagrangeSpace::PdelabWrapper<AluCube3dGridViewType, 1, double, 1>
#endif

#define DISCONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS                                                               \
  Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 1, double, 1>,                      \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<S1dGridPartType, 2, double, 1>,                  \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 1, double, 1>,               \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<Yasp1dGridPartType, 2, double, 1>

#if HAVE_ALUGRID
#define DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM_LOCALFUNCTIONS                                                       \
  Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dGridPartType, 1, double, 1>,             \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluConform2dGridPartType, 2, double, 1>,         \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dGridPartType, 1, double, 1>,         \
      Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<AluSimplex2dGridPartType, 2, double, 1>,         \
      Dune::GDT::DiscontinuousLagrangeSpace::                                                                          \
          FemLocalfunctionsWrapper<AluSimplex3dGridPartType, 1, double, 1> /* <- does not work anymore in 3d */        \
      ,                                                                                                                \
      Dune::GDT::DiscontinuousLagrangeSpace::                                                                          \
          FemLocalfunctionsWrapper<AluSimplex3dGridPartType, 2, double, 1> /* <- does not work anymore in 3d */
#endif // HAVE_ALUGRID

typedef testing::Types<
#if HAVE_DUNE_FEM
    P1_CONTINUOUS_LAGRANGE_SPACES_FEM, Q1_CONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif
#if HAVE_DUNE_FEM_LOCALFUNCTIONS || HAVE_DUNE_PDELAB
    ,
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_FEM_LOCALFUNCTIONS
    P1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS, Q1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM_LOCALFUNCTIONS
#endif
#if HAVE_DUNE_PDELAB
    ,
#endif
#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS
#if HAVE_DUNE_PDELAB
    P1_CONTINUOUS_LAGRANGE_SPACES_PDELAB, Q1_CONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB, Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > P1Q1_Continuous_Lagrange_Spaces;

typedef testing::Types<
#if HAVE_DUNE_FEM
    P1_CONTINUOUS_LAGRANGE_SPACES_FEM, Q1_CONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif
#if HAVE_DUNE_FEM_LOCALFUNCTIONS || HAVE_DUNE_PDELAB
    ,
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_FEM_LOCALFUNCTIONS
    P1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS, Q1_CONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS,
    DISCONTINUOUS_LAGRANGE_SPACES_FEM_LOCALFUNCTIONS
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM_LOCALFUNCTIONS, DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM_LOCALFUNCTIONS
#endif
#if HAVE_DUNE_PDELAB
    ,
#endif
#endif // HAVE_DUNE_FEM_LOCALFUNCTIONS
#if HAVE_DUNE_PDELAB
    P1_CONTINUOUS_LAGRANGE_SPACES_PDELAB, Q1_CONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    P1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB, Q1_CONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > All_Spaces;


TYPED_TEST_CASE(Any_Space, All_Spaces);
TYPED_TEST(Any_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Continuous_Lagrange, P1Q1_Continuous_Lagrange_Spaces);
TYPED_TEST(P1Q1_Continuous_Lagrange, fulfills_interface)
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
