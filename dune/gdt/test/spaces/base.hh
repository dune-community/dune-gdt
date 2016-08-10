// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_TEST_SPACES_BASE_HH
#define DUNE_GDT_TEST_SPACES_BASE_HH

#ifndef DUNE_GDT_TEST_SPACES_RT_CHECK
#define DUNE_GDT_TEST_SPACES_RT_CHECK 0
#endif

#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/basefunctionset/interface.hh>

#include <dune/gdt/test/grids.hh>


// these two should trigger a segfault if copying fails, i.e. the one fixed in 6b3ff6d
template <class Space>
class BaseHolder
{
public:
  BaseHolder(Space s)
    : s_(s)
  {
  }

  const Space& space()
  {
    return s_;
  }

private:
  const Space s_;
};


template <class Space, class Provider>
class DerivedHolder : public BaseHolder<Space>
{
  typedef BaseHolder<Space> BaseType;

public:
  DerivedHolder(Provider& p)
    : BaseType(Space(p.template leaf<Space::part_view_type>()))
  {
  }
};


/**
  * \brief Checks any space derived from SpaceInterface for it's interface compliance, especially concerning CRTP.
  */
template <class SpaceType>
class SpaceBase : public ::testing::Test
{
  typedef typename SpaceType::GridViewType::Grid GridType;
  typedef DSG::Providers::Cube<GridType> ProviderType;

public:
  SpaceBase()
    : grid_provider_(0.0, 1.0, 3u)
    , space_(grid_provider_.template leaf<SpaceType::part_view_type>())
  {
  }

  virtual ~SpaceBase()
  {
  }

  /**
    * \brief Checks the space for it's interface compliance.
    */
  void fulfills_interface() const
  {
    static_assert(Dune::GDT::is_space<SpaceType>::value, "");
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // static checks
    // * as the derived type
    typedef typename SpaceType::Traits Traits;
    typedef typename SpaceType::GridViewType D_GridViewType;
    typedef typename SpaceType::DomainFieldType D_DomainFieldType;
    static const size_t d_dimDomain = SpaceType::dimDomain;
    typedef typename SpaceType::RangeFieldType D_RangeFieldType;
    static const int d_polOrder        = SpaceType::polOrder;
    static const size_t d_dimRange     = SpaceType::dimRange;
    static const size_t d_dimRangeCols = SpaceType::dimRangeCols;
    typedef typename SpaceType::BackendType D_BackendType;
    typedef typename SpaceType::MapperType D_MapperType;
    typedef typename SpaceType::BaseFunctionSetType D_BaseFunctionSetType;
    typedef typename SpaceType::EntityType D_EntityType;
    typedef typename SpaceType::IntersectionType D_IntersectionType;
    typedef typename SpaceType::PatternType D_PatternType;
    typedef typename SpaceType::BoundaryInfoType D_BoundaryInfoType;
    typedef typename SpaceType::CommunicatorType D_CommunicatorType;
    static const bool d_needs_grid_view = SpaceType::needs_grid_view;
    // * as the interface
    typedef SpaceInterface<Traits, d_dimDomain, d_dimRange, d_dimRangeCols> InterfaceType;
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridViewType I_GridViewType;
    typedef typename InterfaceType::DomainFieldType I_DomainFieldType;
    static const size_t i_dimDomain = InterfaceType::dimDomain;
    typedef typename InterfaceType::RangeFieldType I_RangeFieldType;
    static const int i_polOrder        = InterfaceType::polOrder;
    static const size_t i_dimRange     = InterfaceType::dimRange;
    static const size_t i_dimRangeCols = InterfaceType::dimRangeCols;
    typedef typename InterfaceType::BackendType I_BackendType;
    typedef typename InterfaceType::MapperType I_MapperType;
    typedef typename InterfaceType::BaseFunctionSetType I_BaseFunctionSetType;
    typedef typename InterfaceType::EntityType I_EntityType;
    typedef typename InterfaceType::IntersectionType I_IntersectionType;
    typedef typename InterfaceType::PatternType I_PatternType;
    typedef typename InterfaceType::BoundaryInfoType I_BoundaryInfoType;
    typedef typename InterfaceType::CommunicatorType I_CommunicatorType;
    static const bool i_needs_grid_view = InterfaceType::needs_grid_view;
    static_assert(std::is_base_of<InterfaceType, SpaceType>::value, "SpaceType has to be derived from SpaceInterface!");
    static_assert(std::is_same<derived_type, SpaceType>::value, "Types do not match!");
    static_assert(std::is_same<I_GridViewType, D_GridViewType>::value, "Types do not match!");
    static_assert(std::is_same<I_DomainFieldType, D_DomainFieldType>::value, "Types do not match!");
    static_assert(std::is_same<I_RangeFieldType, D_RangeFieldType>::value, "Types do not match!");
    static_assert(std::is_same<I_BackendType, D_BackendType>::value, "Types do not match!");
    static_assert(std::is_same<I_MapperType, D_MapperType>::value, "Types do not match!");
    static_assert(std::is_same<I_BaseFunctionSetType, D_BaseFunctionSetType>::value, "Types do not match!");
    static_assert(std::is_same<I_EntityType, D_EntityType>::value, "Types do not match!");
    static_assert(std::is_same<I_IntersectionType, D_IntersectionType>::value, "Types do not match!");
    static_assert(std::is_same<I_PatternType, D_PatternType>::value, "Types do not match!");
    static_assert(std::is_same<I_BoundaryInfoType, D_BoundaryInfoType>::value, "Types do not match!");
    static_assert(std::is_same<I_CommunicatorType, D_CommunicatorType>::value, "Types do not match!");
    static_assert(std::is_move_constructible<SpaceType>::value, "SpaceType isn't move constructible");
    static_assert(std::is_copy_constructible<SpaceType>::value, "SpaceType isn't copy constructible");
    static_assert(i_dimDomain == d_dimDomain, "Dimensions do not match!");
    static_assert(i_dimRange == d_dimRange, "Dimensions do not match!");
    static_assert(i_dimRangeCols == d_dimRangeCols, "Dimensions do not match!");
    static_assert(i_polOrder == d_polOrder, "Polynomial orders do not match!");
    static_assert(d_needs_grid_view == i_needs_grid_view, "Information do not match!");
    // dynamic checks
    // * as the derived_type
    const D_BackendType& d_backend                 = space_.backend();
    const D_MapperType& d_mapper                   = space_.mapper();
    const D_GridViewType& d_grid_view              = space_.grid_view();
    D_CommunicatorType& d_comm                     = space_.communicator();
    D_PatternType d_pattern                        = space_.compute_pattern();
    D_PatternType d_pattern_view                   = space_.compute_pattern(d_grid_view);
    D_PatternType d_pattern_other                  = space_.compute_pattern(space_);
    D_PatternType d_pattern_view_other             = space_.compute_pattern(d_grid_view, space_);
    D_PatternType d_pattern_volume                 = space_.compute_volume_pattern();
    D_PatternType d_pattern_volume_view            = space_.compute_volume_pattern(d_grid_view);
    D_PatternType d_pattern_volume_other           = space_.compute_volume_pattern(space_);
    D_PatternType d_pattern_volume_view_other      = space_.compute_volume_pattern(d_grid_view, space_);
    D_PatternType d_pattern_face_volume            = space_.compute_face_and_volume_pattern();
    D_PatternType d_pattern_face_volume_view       = space_.compute_face_and_volume_pattern(d_grid_view);
    D_PatternType d_pattern_face_volume_other      = space_.compute_face_and_volume_pattern(space_);
    D_PatternType d_pattern_face_volume_view_other = space_.compute_face_and_volume_pattern(d_grid_view, space_);
    D_PatternType d_pattern_face                   = space_.compute_face_pattern();
    D_PatternType d_pattern_face_view              = space_.compute_face_pattern(d_grid_view);
    D_PatternType d_pattern_face_other             = space_.compute_face_pattern(space_);
    D_PatternType d_pattern_face_view_other        = space_.compute_face_pattern(d_grid_view, space_);
    EXPECT_EQ(d_pattern, d_pattern_other);
    EXPECT_EQ(d_pattern, d_pattern_view);
    EXPECT_EQ(d_pattern, d_pattern_view_other);
    EXPECT_EQ(d_pattern_volume, d_pattern_volume_other);
    EXPECT_EQ(d_pattern_volume, d_pattern_volume_view);
    EXPECT_EQ(d_pattern_volume, d_pattern_volume_view_other);
    EXPECT_EQ(d_pattern_face_volume, d_pattern_face_volume_view);
    EXPECT_EQ(d_pattern_face_volume, d_pattern_face_volume_other);
    EXPECT_EQ(d_pattern_face_volume, d_pattern_face_volume_view_other);
    EXPECT_EQ(d_pattern_face, d_pattern_face_other);
    EXPECT_EQ(d_pattern_face, d_pattern_face_view);
    EXPECT_EQ(d_pattern_face, d_pattern_face_view_other);
    // * as the interface
    const InterfaceType& i_space                   = static_cast<const InterfaceType&>(space_);
    const I_BackendType& i_backend                 = i_space.backend();
    const I_MapperType& i_mapper                   = i_space.mapper();
    const I_GridViewType& i_grid_view              = i_space.grid_view();
    I_CommunicatorType& i_comm                     = i_space.communicator();
    I_PatternType i_pattern                        = i_space.compute_pattern();
    I_PatternType i_pattern_view                   = i_space.compute_pattern(i_grid_view);
    I_PatternType i_pattern_other                  = i_space.compute_pattern(i_space);
    I_PatternType i_pattern_view_other             = i_space.compute_pattern(i_grid_view, i_space);
    I_PatternType i_pattern_volume                 = i_space.compute_volume_pattern();
    I_PatternType i_pattern_volume_view            = i_space.compute_volume_pattern(i_grid_view);
    I_PatternType i_pattern_volume_other           = i_space.compute_volume_pattern(i_space);
    I_PatternType i_pattern_volume_view_other      = i_space.compute_volume_pattern(i_grid_view, i_space);
    I_PatternType i_pattern_face_volume            = i_space.compute_face_and_volume_pattern();
    I_PatternType i_pattern_face_volume_view       = i_space.compute_face_and_volume_pattern(i_grid_view);
    I_PatternType i_pattern_face_volume_other      = i_space.compute_face_and_volume_pattern(i_space);
    I_PatternType i_pattern_face_volume_view_other = i_space.compute_face_and_volume_pattern(i_grid_view, i_space);
    I_PatternType i_pattern_face                   = i_space.compute_face_pattern();
    I_PatternType i_pattern_face_view              = i_space.compute_face_pattern(i_grid_view);
    I_PatternType i_pattern_face_other             = i_space.compute_face_pattern(i_space);
    I_PatternType i_pattern_face_view_other        = i_space.compute_face_pattern(i_grid_view, i_space);
    EXPECT_EQ(&i_backend, &d_backend);
    EXPECT_EQ(&i_mapper, &d_mapper);
    EXPECT_EQ(&i_grid_view, &d_grid_view);
    EXPECT_EQ(&i_comm, &d_comm);
    EXPECT_EQ(i_pattern, d_pattern);
    EXPECT_EQ(i_pattern_other, d_pattern_other);
    EXPECT_EQ(i_pattern_view, d_pattern_view);
    EXPECT_EQ(i_pattern_view_other, d_pattern_view_other);
    EXPECT_EQ(i_pattern_volume, d_pattern_volume);
    EXPECT_EQ(i_pattern_volume_other, d_pattern_volume_other);
    EXPECT_EQ(i_pattern_volume_view, d_pattern_volume_view);
    EXPECT_EQ(i_pattern_volume_view_other, d_pattern_volume_view_other);
    EXPECT_EQ(i_pattern_face_volume, d_pattern_face_volume);
    EXPECT_EQ(i_pattern_face_volume_other, d_pattern_face_volume_other);
    EXPECT_EQ(i_pattern_face_volume_view, d_pattern_face_volume_view);
    EXPECT_EQ(i_pattern_face_volume_view_other, d_pattern_face_volume_view_other);
    EXPECT_EQ(i_pattern_face, d_pattern_face);
    EXPECT_EQ(i_pattern_face_other, d_pattern_face_other);
    EXPECT_EQ(i_pattern_face_view, d_pattern_face_view);
    EXPECT_EQ(i_pattern_face_view_other, d_pattern_face_view_other);
    // walk the grid
    const auto entity_it_end = d_grid_view.template end<0>();
    for (auto entity_it = d_grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const D_EntityType& entity = *entity_it;
      // * as the derived type
      D_BaseFunctionSetType d_base_function_set = space_.base_function_set(entity);
      size_t d_bfs_size                         = d_base_function_set.size();
      EXPECT_EQ(d_bfs_size, d_mapper.numDofs(entity));
      // * as the interface type
      I_BaseFunctionSetType i_base_function_set = i_space.base_function_set(entity);
      size_t i_bfs_size                         = i_base_function_set.size();
      EXPECT_EQ(d_bfs_size, i_bfs_size);
    } // walk the grid
  } // ... fulfills_interface()

  void check_for_correct_copy()
  {
    SpaceType foop(space_);
    auto DUNE_UNUSED(aa) = foop.mapper().size();
    SpaceType cp         = DerivedHolder<SpaceType, ProviderType>(grid_provider_).space();
    auto DUNE_UNUSED(bb) = cp.mapper().size();
  } // ... check_for_correct_copy()

  /**
    * \brief Checks the spaces mapper for it's interface compliance.
    */
  void mapper_fulfills_interface() const
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // static checks
    // * as the derived type
    typedef typename SpaceType::MapperType MapperType;
    typedef typename MapperType::Traits Traits;
    typedef typename MapperType::BackendType D_BackendType;
    // * as the interface
    typedef MapperInterface<Traits> InterfaceType;
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::BackendType I_BackendType;
    static_assert(std::is_base_of<InterfaceType, MapperType>::value,
                  "MapperType has to be derived from MapperInterface!");
    static_assert(std::is_same<derived_type, MapperType>::value, "Types do not match!");
    static_assert(std::is_same<I_BackendType, D_BackendType>::value, "Types do not match!");
    // dynamic checks
    // * as the derived type
    const MapperType& d_mapper     = space_.mapper();
    const D_BackendType& d_backend = d_mapper.backend();
    size_t d_size                  = d_mapper.size();
    size_t d_maxNumDofs            = d_mapper.maxNumDofs();
    // * as the interface type
    const InterfaceType& i_mapper  = static_cast<const InterfaceType&>(d_mapper);
    const D_BackendType& i_backend = i_mapper.backend();
    size_t i_size                  = i_mapper.size();
    size_t i_maxNumDofs            = i_mapper.maxNumDofs();
    EXPECT_EQ(&i_backend, &d_backend);
    EXPECT_EQ(i_size, d_size);
    EXPECT_EQ(i_maxNumDofs, d_maxNumDofs);
    //   walk the grid
    const auto entity_it_end = space_.grid_view().template end<0>();
    for (auto entity_it = space_.grid_view().template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // * as the derived type
      size_t d_numDofs = d_mapper.numDofs(entity);
      DynamicVector<size_t> d_globalIndices(d_numDofs, 0);
      d_mapper.globalIndices(entity, d_globalIndices);
      if (d_globalIndices.size() > d_numDofs)
        DUNE_THROW(Exceptions::index_out_of_range, d_globalIndices.size() << " vs. " << d_numDofs);
      DynamicVector<size_t> d_globalIndices_return = d_mapper.globalIndices(entity);
      EXPECT_EQ(d_globalIndices_return, d_globalIndices);
      // * as the interface
      size_t i_numDofs = i_mapper.numDofs(entity);
      DynamicVector<size_t> i_globalIndices(i_numDofs, 0);
      i_mapper.globalIndices(entity, i_globalIndices);
      DynamicVector<size_t> i_globalIndices_return = i_mapper.globalIndices(entity);
      EXPECT_EQ(i_numDofs, d_numDofs);
      EXPECT_EQ(i_globalIndices, d_globalIndices);
      EXPECT_EQ(i_globalIndices_return, d_globalIndices_return);
      //   walk the local DoFs
      for (size_t ii = 0; ii < d_numDofs; ++ii) {
        // * as the derived type
        size_t d_mapToGlobal = d_mapper.mapToGlobal(entity, ii);
        EXPECT_EQ(d_mapToGlobal, d_globalIndices[ii]);
        // * as the interface
        size_t i_mapToGlobal = i_mapper.mapToGlobal(entity, ii);
        EXPECT_EQ(i_mapToGlobal, d_mapToGlobal);
      } //   walk the local DoFs
    } //   walk the grid
  } // ... mapper_fulfills_interface()

  /**
    * \brief  Checks the spaces basefunctionsets for their interface compliance.
    * \note   We do not check for the functionality enforced by LocalfuntionSetInterface at the moment!
    */
  void basefunctionset_fulfills_interface() const
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    using namespace Dune::GDT;
    // static checks
    // * as the derived type
    typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename BaseFunctionSetType::Traits Traits;
    typedef typename BaseFunctionSetType::BackendType D_BackendType;
    typedef typename BaseFunctionSetType::EntityType D_EntityType;
    typedef typename BaseFunctionSetType::DomainFieldType D_DomainFieldType;
    static const size_t d_dimDomain = BaseFunctionSetType::dimDomain;
    typedef typename BaseFunctionSetType::DomainType D_DomainType;
    typedef typename BaseFunctionSetType::RangeFieldType D_RangeFieldType;
    static const size_t d_dimRange     = BaseFunctionSetType::dimRange;
    static const size_t d_dimRangeCols = BaseFunctionSetType::dimRangeCols;
    typedef typename BaseFunctionSetType::RangeType D_RangeType;
    typedef typename BaseFunctionSetType::JacobianRangeType D_JacobianRangeType;
    static_assert(std::is_same<D_EntityType, typename SpaceType::EntityType>::value, "Types do not match!");
    static_assert(std::is_same<D_DomainFieldType, typename SpaceType::DomainFieldType>::value, "Types do not match!");
    static_assert(std::is_same<D_DomainType, typename SpaceType::DomainType>::value, "Types do not match!");
    static_assert(std::is_same<D_RangeFieldType, typename SpaceType::RangeFieldType>::value, "Types do not match!");
    static_assert(d_dimDomain == SpaceType::dimDomain, "Dimensions do not match!");
    static_assert(d_dimRange == SpaceType::dimRange, "Dimensions do not match!");
    static_assert(d_dimRangeCols == SpaceType::dimRangeCols, "Dimensions do not match!");
    // * as the interface type
    typedef BaseFunctionSetInterface<Traits,
                                     D_DomainFieldType,
                                     d_dimDomain,
                                     D_RangeFieldType,
                                     d_dimRange,
                                     d_dimRangeCols>
        InterfaceType;
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::BackendType I_BackendType;
    typedef typename InterfaceType::EntityType I_EntityType;
    typedef typename InterfaceType::DomainFieldType I_DomainFieldType;
    static const size_t i_dimDomain = InterfaceType::dimDomain;
    typedef typename InterfaceType::DomainType I_DomainType;
    typedef typename InterfaceType::RangeFieldType I_RangeFieldType;
    static const size_t i_dimRange     = InterfaceType::dimRange;
    static const size_t i_dimRangeCols = InterfaceType::dimRangeCols;
    typedef typename InterfaceType::RangeType I_RangeType;
    typedef typename InterfaceType::JacobianRangeType I_JacobianRangeType;
    static_assert(std::is_same<derived_type, BaseFunctionSetType>::value, "Types do not match!");
    static_assert(std::is_same<I_BackendType, D_BackendType>::value, "Types do not match!");
    static_assert(std::is_same<I_EntityType, D_EntityType>::value, "Types do not match!");
    static_assert(std::is_same<I_DomainFieldType, D_DomainFieldType>::value, "Types do not match!");
    static_assert(std::is_same<I_DomainType, D_DomainType>::value, "Types do not match!");
    static_assert(std::is_same<I_RangeFieldType, D_RangeFieldType>::value, "Types do not match!");
    static_assert(std::is_same<I_RangeType, D_RangeType>::value, "Types do not match!");
    static_assert(std::is_same<I_JacobianRangeType, D_JacobianRangeType>::value, "Types do not match!");
    static_assert(i_dimDomain == d_dimDomain, "Dimensions do not match!");
    static_assert(i_dimRange == d_dimRange, "Dimensions do not match!");
    static_assert(i_dimRangeCols == d_dimRangeCols, "Dimensions do not match!");
    // dynamic checks
    // walk the grid
    const auto entity_end_it = space_.grid_view().template end<0>();
    for (auto entity_it = space_.grid_view().template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      // * as the derived type
      BaseFunctionSetType d_base_function_set = space_.base_function_set(entity);
      const D_BackendType& d_backend          = d_base_function_set.backend();
      size_t d_order                          = d_base_function_set.order();
#if DUNE_GDT_TEST_SPACES_RT_CHECK
      EXPECT_GE(d_order, boost::numeric_cast<size_t>(SpaceType::polOrder)); // <- normaly we would expect equality here,
// but the raviart
//    thomas space of order 0 reports order 1 here
#else
      EXPECT_EQ(d_order, boost::numeric_cast<size_t>(SpaceType::polOrder));
#endif
      //   the size has already been checked in fulfills_interface() above
      // * as the interface
      InterfaceType& i_base_function_set = static_cast<InterfaceType&>(d_base_function_set);
      const I_BackendType& i_backend     = i_base_function_set.backend();
      EXPECT_EQ(&d_backend, &i_backend);
      size_t i_order = i_base_function_set.order();
      EXPECT_EQ(i_order, d_order);
    } // walk the grid
  } // ... basefunctionset_fulfills_interface()

protected:
  ProviderType grid_provider_;
  SpaceType space_;
}; // struct SpaceBase

#endif // DUNE_GDT_TEST_SPACES_BASE_HH
