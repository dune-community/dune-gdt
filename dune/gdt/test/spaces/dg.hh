// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)

#ifndef DUNE_GDT_TEST_SPACES_DG_HH
#define DUNE_GDT_TEST_SPACES_DG_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/ranges.hh>

#include "base.hh"


template <class SpaceType>
class DG_Space : public SpaceBase<SpaceType>
{
};


template <class SpaceType>
struct P1Q1_DG_Space : public SpaceBase<SpaceType>
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  static const size_t dimDomain = GridType::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  static std::vector<DomainFieldType> convert_vector(const DomainType& source)
  {
    std::vector<DomainFieldType> ret(dimDomain, DomainFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii)
      ret[ii]      = source[ii];
    return ret;
  }

  void maps_correctly()
  {
    using namespace Dune::Stuff;
    // walk the grid to create a map of all vertices
    std::map<std::vector<DomainFieldType>, std::pair<std::set<size_t>, size_t>> vertex_to_indices_map;
    const auto entity_end_it = this->space_.grid_view().template end<0>();
    for (auto entity_it = this->space_.grid_view().template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      for (auto cc : DSC::valueRange(entity.subEntities(dimDomain))) {
        const auto vertex                                          = entity.template subEntity<dimDomain>(cc);
        const DomainType vertex_center                             = vertex.geometry().center();
        vertex_to_indices_map[convert_vector(vertex_center)].first = std::set<size_t>();
        ++vertex_to_indices_map[convert_vector(vertex_center)].second;
      }
    }
    // walk the grid again to find all DoF ids
    for (auto entity_it = this->space_.grid_view().template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity        = *entity_it;
      const size_t num_vertices = entity.subEntities(dimDomain);
      const auto basis          = this->space_.base_function_set(entity);
      EXPECT_EQ(basis.size(), num_vertices);
      for (size_t cc = 0; cc < num_vertices; ++cc) {
        const auto vertex              = entity.template subEntity<dimDomain>(boost::numeric_cast<int>(cc));
        const DomainType vertex_center = vertex.geometry().center();
        // find the local basis function which corresponds to this vertex
        const auto basis_values = basis.evaluate(entity.geometry().local(vertex_center));
        EXPECT_EQ(basis_values.size(), num_vertices);
        size_t ones            = 0;
        size_t zeros           = 0;
        size_t failures        = 0;
        size_t local_DoF_index = 0;
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          if (Common::FloatCmp::eq(basis_values[ii][0], typename SpaceType::RangeFieldType(1))) {
            local_DoF_index = ii;
            ++ones;
          } else if (Common::FloatCmp::eq(basis_values[ii][0] + 1, typename SpaceType::RangeFieldType(1)))
            ++zeros;
          else
            ++failures;
        }
        if (ones != 1 || zeros != (basis.size() - 1) || failures > 0) {
          std::stringstream ss;
          ss << "ones = " << ones << ", zeros = " << zeros << ", failures = " << failures
             << ", num_vertices = " << num_vertices << ", entity " << this->space_.grid_view().indexSet().index(entity)
             << ", vertex " << cc << ": [ " << vertex_center << "], ";
          Common::print(basis_values, "basis_values", ss);
          EXPECT_TRUE(false) << ss.str();
        }
        // now we know that the local DoF index of this vertex is ii
        const size_t global_DoF_index = this->space_.mapper().mapToGlobal(entity, local_DoF_index);
        vertex_to_indices_map[convert_vector(vertex_center)].first.insert(global_DoF_index);
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
      EXPECT_EQ(global_DoF_id, count);
      ++count;
    }
    for (const auto& entry : vertex_to_indices_map) {
      const auto vertex_ids               = entry.second.first;
      size_t number_of_associated_DoF_ids = vertex_ids.size();
      size_t number_of_adjacent_entitys   = entry.second.second;
      EXPECT_EQ(number_of_associated_DoF_ids, number_of_adjacent_entitys);
    }
  } // ... maps_correctly()
}; // struct P1Q1_DG_Space


#endif // DUNE_GDT_TEST_SPACES_DG_HH
