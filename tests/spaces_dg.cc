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

template <class SpaceType>
struct P1Q1_Discontinuous_Lagrange : public ::testing::Test, public ::SpaceTestBase<SpaceType>
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
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
    // check that each vertex has the appropiate number of local DoF ids and that the numbering is consecutive
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids = entry.second;
      for (auto vertex_ids_it = vertex_ids.begin(); vertex_ids_it != vertex_ids.end(); ++vertex_ids_it)
        global_DoF_indices.insert(*vertex_ids_it);
    }
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      if (global_DoF_id != count)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, "count = " << count << ", global_DoF_id = " << global_DoF_id);
      ++count;
    }
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_coordinates = entry.first;
      const auto vertex_ids         = entry.second;
      size_t vertex_boundary_count = 0;
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        if (Common::FloatCmp::eq(vertex_coordinates[ii], DomainFieldType(1))
            || Common::FloatCmp::eq(vertex_coordinates[ii], DomainFieldType(0)))
          ++vertex_boundary_count;
      }
      /* we are on a cubic grid, so each vertex in the interior of the grid should have pow(2, dimDomain) adjacent
      entitys. If the vertex is on a face of the grid (vertex_boundary_count = 1), the number halves, if it is on the
      intersection of exactly two faces (vertex_boundary_count = 2), it halves again, ... , if the vertex is also an
      vertex of the grid domain (vertex_boundary_count = dimDomain), it has only 1 adjacent entity */
      size_t adjacent_entitys_expected = pow(2, dimDomain - vertex_boundary_count);
      size_t adjacent_entitys = vertex_ids.size();
      if (adjacent_entitys != adjacent_entitys_expected)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error,
                              "Vertex has only " << adjacent_entitys << "adjacent entitys, should have "
                                                 << adjacent_entitys_expected);
    }

  } // ... maps_correctly()
}; // struct P1Q1_Discontinuous_Lagrange


#if HAVE_DUNE_FEM
#define Qk_DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 1, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 2, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 2, double, 1>

#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 1, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S2dLeafGridPartType, 1, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S3dLeafGridPartType, 1, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp2dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp3dLeafGridPartType, 1, double, 1>

#if HAVE_ALUGRID
#define Pk_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 2, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 2, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 2, double, 1>

#define P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 1, double, 1>

#define P1D_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                  \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 1, double, 2>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 1, double, 2>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 1, double, 3>

#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube3dLeafGridPartType, 1, double, 1>
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_PDELAB
#define Qk_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 1, double, 1>,                            \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 2, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 1, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 2, double, 1>


#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 1, double, 1>,                            \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S2dLeafGridViewType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S3dLeafGridViewType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 1, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp2dLeafGridViewType, 1, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp3dLeafGridViewType, 1, double, 1>

#if HAVE_ALUGRID
#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube3dLeafGridViewType, 1, double, 1>
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB

typedef testing::Types<
#if 0 // HAVE_DUNE_FEM
                        P1_DISCONTINUOUS_LAGRANGE_SPACES_FEM
                      , Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM
                      , Qk_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
                      , P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
                      , Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
                      , Pk_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
                      , P1D_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif
#if HAVE_DUNE_PDELAB
                      ,
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB, Qk_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > All_Spaces;


typedef testing::Types<
#if 0 // HAVE_DUNE_FEM
                        Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
                        Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif
#if HAVE_DUNE_PDELAB
                      ,
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > Q1_Spaces;

typedef testing::Types<
#if HAVE_DUNE_FEM && HAVE_ALUGRID
    P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, P1D_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif // HAVE_DUNE_FEM && HAVE_ALUGRID
    > P1_Spaces;


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

TYPED_TEST_CASE(P1Q1_Discontinuous_Lagrange, Q1_Spaces);
TYPED_TEST(P1Q1_Discontinuous_Lagrange, maps_correctly)
{
  this->maps_correctly();
}

// TYPED_TEST_CASE(P1Q1_Discontinuous_Lagrange, P1_Spaces);
// TYPED_TEST(P1Q1_Discontinuous_Lagrange, maps_correctly)
//{
//  this->maps_correctly();
//}


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
