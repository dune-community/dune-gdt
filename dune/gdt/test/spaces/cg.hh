// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_TEST_SPACES_CG_HH
#define DUNE_GDT_TEST_SPACES_CG_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/print.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/cg/dune-fem-wrapper.hh>
#include <dune/gdt/spaces/cg/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/basefunctionset/interface.hh>

#include "base.hh"


template <class SpaceType>
class CG_Space : public SpaceBase<SpaceType>
{
};


template <class SpaceType>
struct P1Q1_CG_Space : public SpaceBase<SpaceType>
{
  typedef typename SpaceType::GridLayerType GridLayerType;
  using GridType = Dune::XT::Grid::extract_grid_t<typename SpaceType::GridLayerType>;
  static const size_t dimDomain = GridType::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  static std::vector<DomainFieldType> convert_vector(const DomainType& source)
  {
    std::vector<DomainFieldType> ret(dimDomain, DomainFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret[ii] = source[ii];
    return ret;
  }

  template <class T, size_t d, size_t r, size_t rC>
  void matches_signature(const Dune::GDT::CgSpaceInterface<T, d, r, rC>& /*space*/)
  {
    static_assert(std::is_same<typename SpaceType::Traits, T>::value, "");
    static_assert(d == SpaceType::dimDomain, "");
    static_assert(r == SpaceType::dimRange, "");
    static_assert(rC == SpaceType::dimRangeCols, "");
  }

  void fulfills_continuous_interface()
  {

    matches_signature(this->space_);
    const auto entity_ptr = this->space_.grid_layer().template begin<0>();
    const auto& entity = *entity_ptr;
    const auto basis = this->space_.base_function_set(entity);
    std::vector<DomainType> lagrange_points = this->space_.lagrange_points(entity);
    EXPECT_EQ(lagrange_points.size(), basis.size());
    typedef Dune::XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType> IntersectionType;
    typedef typename SpaceType::RangeFieldType RangeFieldType DUNE_UNUSED;
    Dune::XT::Grid::AllDirichletBoundaryInfo<IntersectionType> boundary_info;
    std::set<size_t> local_dirichlet_DoFs DUNE_UNUSED = this->space_.local_dirichlet_DoFs(entity, boundary_info);
    Dune::GDT::DirichletConstraints<IntersectionType> dirichlet_constraints_set(
        boundary_info, this->space_.mapper().size(), true);
    Dune::GDT::DirichletConstraints<IntersectionType> dirichlet_constraints_clear(
        boundary_info, this->space_.mapper().size(), false);
    this->space_.local_constraints(entity, dirichlet_constraints_set);
    this->space_.local_constraints(entity, dirichlet_constraints_clear);
  }

  void maps_correctly()
  {
    using namespace Dune::XT;
    // walk the grid to create a map of all vertices
    std::map<std::vector<DomainFieldType>, std::set<size_t>> vertex_to_indices_map;
    const auto entity_end_it = this->space_.grid_layer().template end<0>();
    for (auto entity_it = this->space_.grid_layer().template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      for (auto cc : Dune::XT::Common::value_range(entity.subEntities(dimDomain))) {
        const auto vertex = entity.template subEntity<dimDomain>(cc);
        const DomainType vertex_center = vertex.geometry().center();
        vertex_to_indices_map[convert_vector(vertex_center)] = std::set<size_t>();
      }
    }

    // walk the grid again to find all DoF ids
    auto functor = [&](const Dune::XT::Grid::extract_entity_t<GridLayerType>& entity) {
      const size_t num_vertices = entity.subEntities(dimDomain);
      const auto basis = this->space_.base_function_set(entity);
      EXPECT_EQ(basis.size(), num_vertices);
      for (size_t cc = 0; cc < num_vertices; ++cc) {
        const auto vertex = entity.template subEntity<dimDomain>(boost::numeric_cast<int>(cc));
        const DomainType vertex_center = vertex.geometry().center();
        // find the local basis function which corresponds to this vertex
        const auto basis_values = basis.evaluate(entity.geometry().local(vertex_center));
        EXPECT_EQ(basis_values.size(), num_vertices);
        size_t ones = 0;
        size_t zeros = 0;
        size_t failures = 0;
        size_t local_DoF_index = 0;
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          if (Common::FloatCmp::eq(basis_values[ii][0], typename SpaceType::RangeFieldType(1))) {
            local_DoF_index = ii;
            ++ones;
          } else if (Common::FloatCmp::eq(basis_values[ii][0] + 1.0, typename SpaceType::RangeFieldType(1)))
            ++zeros;
          else
            ++failures;
        }
        if (ones != 1 || zeros != (basis.size() - 1) || failures > 0) {
          std::stringstream ss;
          ss << "ones = " << ones << ", zeros = " << zeros << ", failures = " << failures
             << ", num_vertices = " << num_vertices << ", entity " << this->space_.grid_layer().indexSet().index(entity)
             << ", vertex " << cc << ": [ " << vertex_center << "], ";
          Dune::XT::Common::print(basis_values, "basis_values", ss);
          EXPECT_TRUE(false) << ss.str();
        }
        // now we know that the local DoF index of this vertex is ii
        const size_t global_DoF_index = this->space_.mapper().mapToGlobal(entity, local_DoF_index);
        vertex_to_indices_map[convert_vector(vertex_center)].insert(global_DoF_index);
      }
    };
    Dune::XT::Grid::Walker<GridLayerType> walker(this->space_.grid_layer());
    walker.append(functor);
    walker.walk();
    // check that all vertices have indeed one and only one global DoF id and that the numbering is consecutive
    std::set<size_t> global_DoF_indices;
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids = entry.second;
      EXPECT_EQ(vertex_ids.size(), 1);
      global_DoF_indices.insert(*(vertex_ids.begin()));
    }
    EXPECT_EQ(vertex_to_indices_map.size(), global_DoF_indices.size());
    size_t count = 0;
    for (const auto& global_DoF_id : global_DoF_indices) {
      EXPECT_EQ(global_DoF_id, count);
      ++count;
    }
  } // ... maps_correctly()
}; // struct P1Q1_CG_Space


#endif // DUNE_GDT_TEST_SPACES_CG_HH
