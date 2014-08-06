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
#include <utility>

#include <dune/common/typetraits.hh>
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/print.hh>

#include <dune/gdt/playground/spaces/discontinuouslagrange/fem.hh>
#include <dune/gdt/playground/spaces/discontinuouslagrange/pdelab.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


template <class SpaceType>
class P1Q1_Space : public ::testing::Test, public SpaceTestBase<SpaceType>
{
};

template <class SpaceType>
class P2Q2_Space : public ::testing::Test, public SpaceTestBase<SpaceType>
{
};

template <class SpaceType>
class P3Q3_Space : public ::testing::Test, public SpaceTestBase<SpaceType>
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
    auto& grid_ptr            = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid_ptr);
    const SpaceType space(grid_part_view);
    // walk the grid to create a map of all vertices
    std::map<std::vector<DomainFieldType>, std::pair<std::set<size_t>, size_t>> vertex_to_indices_map;
    const auto entity_end_it = grid_part_view->template end<0>();
    for (auto entity_it = grid_part_view->template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      for (int cc = 0; cc < entity.template count<dimDomain>(); ++cc) {
        const auto vertex_ptr   = entity.template subEntity<dimDomain>(cc);
        const DomainType vertex = vertex_ptr->geometry().center();
        vertex_to_indices_map[convert_vector(vertex)].first = std::set<size_t>();
        ++vertex_to_indices_map[convert_vector(vertex)].second;
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
        vertex_to_indices_map[convert_vector(vertex)].first.insert(global_DoF_index);
      }
    }
    // check that each vertex has the appropiate number of associated DoF ids and that the numbering is consecutive
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids = entry.second.first;
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
      const auto vertex_coordinates       = entry.first;
      const auto vertex_ids               = entry.second.first;
      size_t number_of_associated_DoF_ids = vertex_ids.size();
      size_t number_of_adjacent_entitys = entry.second.second;
      if (number_of_associated_DoF_ids != number_of_adjacent_entitys)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error,
                              "Vertex has only " << number_of_associated_DoF_ids << "associated DoF_ids, should have "
                                                 << number_of_adjacent_entitys);
    }

  } // ... maps_correctly()
}; // struct P1Q1_Discontinuous_Lagrange


#if HAVE_DUNE_FEM
#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 1, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S2dLeafGridPartType, 1, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S3dLeafGridPartType, 1, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp2dLeafGridPartType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp3dLeafGridPartType, 1, double, 1>

#define Q2_DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 2, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S2dLeafGridPartType, 2, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S3dLeafGridPartType, 2, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 2, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp2dLeafGridPartType, 2, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp3dLeafGridPartType, 2, double, 1>

#define Q3_DISCONTINUOUS_LAGRANGE_SPACES_FEM                                                                           \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S1dLeafGridPartType, 3, double, 1>,                               \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S2dLeafGridPartType, 3, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<S3dLeafGridPartType, 3, double, 1>,                           \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp1dLeafGridPartType, 3, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp2dLeafGridPartType, 3, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<Yasp3dLeafGridPartType, 3, double, 1>

#if HAVE_ALUGRID
#define P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 1, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 1, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 1, double, 1>

#define P2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 2, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 2, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 2, double, 1>

#define P3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluConform2dLeafGridPartType, 3, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex2dLeafGridPartType, 3, double, 1>,                  \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluSimplex3dLeafGridPartType, 3, double, 1>

#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube2dLeafGridPartType, 1, double, 1>,                         \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube3dLeafGridPartType, 1, double, 1>

#define Q2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube2dLeafGridPartType, 2, double, 1>,                         \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube3dLeafGridPartType, 2, double, 1>

#define Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM                                                                   \
  Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube2dLeafGridPartType, 3, double, 1>,                         \
      Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<AluCube3dLeafGridPartType, 3, double, 1>
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_PDELAB
#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 1, double, 1>,                            \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S2dLeafGridViewType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S3dLeafGridViewType, 1, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 1, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp2dLeafGridViewType, 1, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp3dLeafGridViewType, 1, double, 1>

#define Q2_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 2, double, 1>,                            \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S2dLeafGridViewType, 2, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S3dLeafGridViewType, 2, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 2, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp2dLeafGridViewType, 2, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp3dLeafGridViewType, 2, double, 1>

#define Q3_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB                                                                        \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S1dLeafGridViewType, 3, double, 1>,                            \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S2dLeafGridViewType, 3, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<S3dLeafGridViewType, 3, double, 1>,                        \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp1dLeafGridViewType, 3, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp2dLeafGridViewType, 3, double, 1>,                     \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<Yasp3dLeafGridViewType, 3, double, 1>

#if HAVE_ALUGRID
#define Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube2dLeafGridViewType, 1, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube3dLeafGridViewType, 1, double, 1>

#define Q2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube2dLeafGridViewType, 2, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube3dLeafGridViewType, 2, double, 1>

#define Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB                                                                \
  Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube2dLeafGridViewType, 3, double, 1>,                      \
      Dune::GDT::Spaces::DiscontinuousLagrange::PdelabBased<AluCube3dLeafGridViewType, 3, double, 1>
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_PDELAB


typedef testing::Types<
#if HAVE_DUNE_FEM
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q1_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > P1Q1_Spaces;

typedef testing::Types<
#if HAVE_DUNE_FEM
    Q2_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    ,
    Q2_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q2_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > P2Q2_Spaces;

typedef testing::Types<
#if HAVE_DUNE_FEM
    Q3_DISCONTINUOUS_LAGRANGE_SPACES_FEM
#if HAVE_ALUGRID
    ,
    P3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM, Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_FEM
#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
    ,
    Q3_DISCONTINUOUS_LAGRANGE_SPACES_PDELAB
#if HAVE_ALUGRID
    ,
    Q3_DISCONTINUOUS_LAGRANGE_SPACES_ALUGRID_PDELAB
#endif
#endif // HAVE_DUNE_PDELAB
    > P3Q3_Spaces;


TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P2Q2_Space, P2Q2_Spaces);
TYPED_TEST(P2Q2_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, fulfills_interface)
{
  this->fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(P2Q2_Space, P2Q2_Spaces);
TYPED_TEST(P2Q2_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Space, P1Q1_Spaces);
TYPED_TEST(P1Q1_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}

TYPED_TEST_CASE(P2Q2_Space, P2Q2_Spaces);
TYPED_TEST(P2Q2_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}

TYPED_TEST_CASE(P3Q3_Space, P3Q3_Spaces);
TYPED_TEST(P3Q3_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}

TYPED_TEST_CASE(P1Q1_Discontinuous_Lagrange, P1Q1_Spaces);
TYPED_TEST(P1Q1_Discontinuous_Lagrange, maps_correctly)
{
  this->maps_correctly();
}

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
