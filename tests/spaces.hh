// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <type_traits>

#if HAVE_ALUGRID
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/alugrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/basefunctionset/interface.hh>


#define SGRID_TYPES(dim)                                                                                               \
  typedef Dune::SGrid<dim, dim> S##dim##dGridType;                                                                     \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S##dim##dGridType, false>::Type S##dim##dLeafGridPartType;  \
  typedef                                                                                                              \
      typename Dune::GDT::SpaceTools::LevelGridPartView<S##dim##dGridType, false>::Type S##dim##dLevelGridPartType;    \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView<S##dim##dGridType, true>::Type S##dim##dLeafGridViewType;   \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView<S##dim##dGridType, true>::Type S##dim##dLevelGridViewType;
SGRID_TYPES(1)
SGRID_TYPES(2)
SGRID_TYPES(3)
#undef SGRID_TYPES

#define YASPGRID_TYPES(dim)                                                                                            \
  typedef Dune::YaspGrid<dim> Yasp##dim##dGridType;                                                                    \
  typedef typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp##dim##dGridType, false>::Type                          \
      Yasp##dim##dLeafGridPartType;                                                                                    \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView<Yasp##dim##dGridType, false>::Type                         \
      Yasp##dim##dLevelGridPartType;                                                                                   \
  typedef                                                                                                              \
      typename Dune::GDT::SpaceTools::LeafGridPartView<Yasp##dim##dGridType, true>::Type Yasp##dim##dLeafGridViewType; \
  typedef typename Dune::GDT::SpaceTools::LevelGridPartView<Yasp##dim##dGridType, true>::Type                          \
      Yasp##dim##dLevelGridViewType;
YASPGRID_TYPES(1)
YASPGRID_TYPES(2)
YASPGRID_TYPES(3)
#undef YASPGRID_TYPES


#if HAVE_ALUGRID
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> AluConform2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, false>::Type AluConform2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, true>::Type AluConform2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, true>::Type AluConform2dLevelGridViewType;
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> AluSimplex2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLevelGridViewType;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming> AluSimplex3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLevelGridViewType;
typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> AluCube2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, false>::Type AluCube2dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, false>::Type AluCube2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube2dGridType, true>::Type AluCube2dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube2dGridType, true>::Type AluCube2dLevelGridViewType;
typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming> AluCube3dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, false>::Type AluCube3dLeafGridPartType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, false>::Type AluCube3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluCube3dGridType, true>::Type AluCube3dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LevelGridPartView<AluCube3dGridType, true>::Type AluCube3dLevelGridViewType;
#endif


/**
  * \brief Checks any space derived from SpaceInterface for it's interface compliance, especially concerning CRTP.
  */
template <class SpaceType>
class SpaceTestBase
{
  typedef typename SpaceType::GridViewType::Grid GridType;
  typedef DSG::Providers::Cube<GridType> ProviderType;

public:
  ~SpaceTestBase()
  {
  }

  SpaceTestBase()
    : grid_provider_(0.0, 1.0, 3u)
  {
    using namespace Dune;
    using namespace GDT;
    const auto grid_part_view = SpaceTools::GridPartView<SpaceType>::create_leaf(grid_provider_.grid());
    space_                    = std::unique_ptr<SpaceType>(new SpaceType(grid_part_view));
  }

  /**
    * \brief Checks the space for it's interface compliance.
    */
  void fulfills_interface() const
  {
    using namespace Dune;
    using namespace GDT;
    using namespace Stuff;
    if (!space_)
      DUNE_THROW_COLORFULLY(Exceptions::internal_error, "");
    // static checks
    // * as the derived type
    typedef typename SpaceType::Traits Traits;
    typedef typename SpaceType::GridViewType D_GridViewType;
    typedef typename SpaceType::DomainFieldType D_DomainFieldType;
    static const unsigned int d_dimDomain = SpaceType::dimDomain;
    typedef typename SpaceType::RangeFieldType D_RangeFieldType;
    static const unsigned int d_dimRange     = SpaceType::dimRange;
    static const unsigned int d_dimRangeCols = SpaceType::dimRangeCols;
    static const int d_polOrder              = SpaceType::polOrder;
    typedef typename SpaceType::BackendType D_BackendType;
    typedef typename SpaceType::MapperType D_MapperType;
    typedef typename SpaceType::BaseFunctionSetType D_BaseFunctionSetType;
    typedef typename SpaceType::EntityType D_EntityType;
    typedef typename SpaceType::IntersectionType D_IntersectionType;
    typedef typename SpaceType::PatternType D_PatternType;
    typedef typename SpaceType::BoundaryInfoType D_BoundaryInfoType;
    static const bool d_needs_grid_view = SpaceType::needs_grid_view;
    // * as the interface
    typedef SpaceInterface<Traits> InterfaceType;
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::GridViewType I_GridViewType;
    typedef typename InterfaceType::DomainFieldType I_DomainFieldType;
    static const unsigned int i_dimDomain = InterfaceType::dimDomain;
    typedef typename InterfaceType::RangeFieldType I_RangeFieldType;
    static const unsigned int i_dimRange     = InterfaceType::dimRange;
    static const unsigned int i_dimRangeCols = InterfaceType::dimRangeCols;
    static const int i_polOrder              = InterfaceType::polOrder;
    typedef typename InterfaceType::BackendType I_BackendType;
    typedef typename InterfaceType::MapperType I_MapperType;
    typedef typename InterfaceType::BaseFunctionSetType I_BaseFunctionSetType;
    typedef typename InterfaceType::EntityType I_EntityType;
    typedef typename InterfaceType::IntersectionType I_IntersectionType;
    typedef typename InterfaceType::PatternType I_PatternType;
    typedef typename InterfaceType::BoundaryInfoType I_BoundaryInfoType;
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
    static_assert(i_dimDomain == d_dimDomain, "Dimensions do not match!");
    static_assert(i_dimRange == d_dimRange, "Dimensions do not match!");
    static_assert(i_dimRangeCols == d_dimRangeCols, "Dimensions do not match!");
    static_assert(i_polOrder == d_polOrder, "Polynomial orders do not match!");
    static_assert(d_needs_grid_view == i_needs_grid_view, "Information do not match!");
    // dynamic checks
    // * as the derived_type
    const D_BackendType& d_backend                           = space_->backend();
    const D_MapperType& d_mapper                             = space_->mapper();
    const std::shared_ptr<const D_GridViewType>& d_grid_view = space_->grid_view();
    D_PatternType d_pattern                                  = space_->compute_pattern();
    D_PatternType d_pattern_view                             = space_->compute_pattern(*d_grid_view);
    D_PatternType d_pattern_other                            = space_->compute_pattern(*space_);
    D_PatternType d_pattern_view_other                       = space_->compute_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_volume                           = space_->compute_volume_pattern();
    D_PatternType d_pattern_volume_view                      = space_->compute_volume_pattern(*d_grid_view);
    D_PatternType d_pattern_volume_other                     = space_->compute_volume_pattern(*space_);
    D_PatternType d_pattern_volume_view_other                = space_->compute_volume_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_face_volume                      = space_->compute_face_and_volume_pattern();
    D_PatternType d_pattern_face_volume_view                 = space_->compute_face_and_volume_pattern(*d_grid_view);
    D_PatternType d_pattern_face_volume_other                = space_->compute_face_and_volume_pattern(*space_);
    D_PatternType d_pattern_face_volume_view_other           = space_->compute_face_and_volume_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_face                             = space_->compute_face_pattern();
    D_PatternType d_pattern_face_view                        = space_->compute_face_pattern(*d_grid_view);
    D_PatternType d_pattern_face_other                       = space_->compute_face_pattern(*space_);
    D_PatternType d_pattern_face_view_other = space_->compute_face_pattern(*d_grid_view, *space_);
    if (d_pattern != d_pattern_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern != d_pattern_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern != d_pattern_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_volume != d_pattern_volume_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_volume != d_pattern_volume_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_volume != d_pattern_volume_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face_volume != d_pattern_face_volume_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face_volume != d_pattern_face_volume_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face_volume != d_pattern_face_volume_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face != d_pattern_face_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face != d_pattern_face_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (d_pattern_face != d_pattern_face_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    // * as the interface
    const InterfaceType& i_space                             = static_cast<const InterfaceType&>(*space_);
    const I_BackendType& i_backend                           = i_space.backend();
    const D_MapperType& i_mapper                             = i_space.mapper();
    const std::shared_ptr<const I_GridViewType>& i_grid_view = i_space.grid_view();
    I_PatternType i_pattern                                  = i_space.compute_pattern();
    I_PatternType i_pattern_view                             = i_space.compute_pattern(*i_grid_view);
    I_PatternType i_pattern_other                            = i_space.compute_pattern(i_space);
    I_PatternType i_pattern_view_other                       = i_space.compute_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_volume                           = i_space.compute_volume_pattern();
    I_PatternType i_pattern_volume_view                      = i_space.compute_volume_pattern(*i_grid_view);
    I_PatternType i_pattern_volume_other                     = i_space.compute_volume_pattern(i_space);
    I_PatternType i_pattern_volume_view_other                = i_space.compute_volume_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_face_volume                      = i_space.compute_face_and_volume_pattern();
    I_PatternType i_pattern_face_volume_view                 = i_space.compute_face_and_volume_pattern(*i_grid_view);
    I_PatternType i_pattern_face_volume_other                = i_space.compute_face_and_volume_pattern(i_space);
    I_PatternType i_pattern_face_volume_view_other           = i_space.compute_face_and_volume_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_face                             = i_space.compute_face_pattern();
    I_PatternType i_pattern_face_view                        = i_space.compute_face_pattern(*i_grid_view);
    I_PatternType i_pattern_face_other                       = i_space.compute_face_pattern(i_space);
    I_PatternType i_pattern_face_view_other = i_space.compute_face_pattern(*i_grid_view, i_space);
    if (&i_backend != &d_backend)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (&i_mapper != &d_mapper)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (&i_grid_view != &d_grid_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern != d_pattern)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_other != d_pattern_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_view != d_pattern_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_view_other != d_pattern_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_volume != d_pattern_volume)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_volume_other != d_pattern_volume_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_volume_view != d_pattern_volume_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_volume_view_other != d_pattern_volume_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_volume != d_pattern_face_volume)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_volume_other != d_pattern_face_volume_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_volume_view != d_pattern_face_volume_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_volume_view_other != d_pattern_face_volume_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face != d_pattern_face)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_other != d_pattern_face_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_view != d_pattern_face_view)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_pattern_face_view_other != d_pattern_face_view_other)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    // walk the grid
    const auto entity_it_end = d_grid_view->template end<0>();
    for (auto entity_it = d_grid_view->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const D_EntityType& entity = *entity_it;
      // * s the derived type
      D_BaseFunctionSetType d_base_function_set = space_->base_function_set(entity);
      size_t d_bfs_size = d_base_function_set.size();
      if (d_bfs_size != d_mapper.numDofs(entity))
        DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range, d_bfs_size << " vs. " << d_mapper.numDofs(entity));
      I_BaseFunctionSetType i_base_function_set = i_space.base_function_set(entity);
      size_t i_bfs_size = i_base_function_set.size();
      if (d_bfs_size != i_bfs_size)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    } // walk the grid
  } // ... fulfills_interface()

  /**
    * \brief Checks the spaces mapper for it's interface compliance.
    */
  void mapper_fulfills_interface() const
  {
    using namespace Dune;
    using namespace GDT;
    using namespace Stuff;
    if (!space_)
      DUNE_THROW_COLORFULLY(Exceptions::internal_error, "");
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
    const MapperType& d_mapper     = space_->mapper();
    const D_BackendType& d_backend = d_mapper.backend();
    size_t d_size                  = d_mapper.size();
    size_t d_maxNumDofs            = d_mapper.maxNumDofs();
    // * as the interface type
    const InterfaceType& i_mapper  = static_cast<const InterfaceType&>(d_mapper);
    const D_BackendType& i_backend = i_mapper.backend();
    size_t i_size                  = i_mapper.size();
    size_t i_maxNumDofs = i_mapper.maxNumDofs();
    if (&i_backend != &d_backend)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_size != d_size)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    if (i_maxNumDofs != d_maxNumDofs)
      DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    //   walk the grid
    const auto entity_it_end = space_->grid_view()->template end<0>();
    for (auto entity_it = space_->grid_view()->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // * as the derived type
      size_t d_numDofs = d_mapper.numDofs(entity);
      DynamicVector<size_t> d_globalIndices(d_numDofs, 0);
      d_mapper.globalIndices(entity, d_globalIndices);
      if (d_globalIndices.size() > d_numDofs)
        DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range, d_globalIndices.size() << " vs. " << d_numDofs);
      DynamicVector<size_t> d_globalIndices_return = d_mapper.globalIndices(entity);
      if (d_globalIndices_return != d_globalIndices)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
      // * as the interface
      size_t i_numDofs = i_mapper.numDofs(entity);
      DynamicVector<size_t> i_globalIndices(i_numDofs, 0);
      i_mapper.globalIndices(entity, i_globalIndices);
      DynamicVector<size_t> i_globalIndices_return = i_mapper.globalIndices(entity);
      if (i_numDofs != d_numDofs)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
      if (i_globalIndices != d_globalIndices)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
      if (i_globalIndices_return != d_globalIndices_return)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
      //   walk the local DoFs
      for (size_t ii = 0; ii < d_numDofs; ++ii) {
        // * as the derived type
        size_t d_mapToGlobal = d_mapper.mapToGlobal(entity, ii);
        if (d_mapToGlobal != d_globalIndices[ii])
          DUNE_THROW_COLORFULLY(Exceptions::index_out_of_range, d_mapToGlobal << " vs. " << d_globalIndices[ii]);
        // * as the interface
        size_t i_mapToGlobal = i_mapper.mapToGlobal(entity, ii);
        if (i_mapToGlobal != d_mapToGlobal)
          DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
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
    using namespace GDT;
    using namespace Stuff;
    if (!space_)
      DUNE_THROW_COLORFULLY(Exceptions::internal_error, "");
    // static checks
    // * as the derived type
    typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename BaseFunctionSetType::Traits Traits;
    typedef typename BaseFunctionSetType::BackendType D_BackendType;
    typedef typename BaseFunctionSetType::EntityType D_EntityType;
    typedef typename BaseFunctionSetType::DomainFieldType D_DomainFieldType;
    static const unsigned int d_dimDomain = BaseFunctionSetType::dimDomain;
    typedef typename BaseFunctionSetType::DomainType D_DomainType;
    typedef typename BaseFunctionSetType::RangeFieldType D_RangeFieldType;
    static const unsigned int d_dimRange     = BaseFunctionSetType::dimRange;
    static const unsigned int d_dimRangeCols = BaseFunctionSetType::dimRangeCols;
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
                                     d_dimRangeCols> InterfaceType;
    typedef typename InterfaceType::derived_type derived_type;
    typedef typename InterfaceType::BackendType I_BackendType;
    typedef typename InterfaceType::EntityType I_EntityType;
    typedef typename InterfaceType::DomainFieldType I_DomainFieldType;
    static const unsigned int i_dimDomain = InterfaceType::dimDomain;
    typedef typename InterfaceType::DomainType I_DomainType;
    typedef typename InterfaceType::RangeFieldType I_RangeFieldType;
    static const unsigned int i_dimRange     = InterfaceType::dimRange;
    static const unsigned int i_dimRangeCols = InterfaceType::dimRangeCols;
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
    const auto entity_end_it = space_->grid_view()->template end<0>();
    for (auto entity_it = space_->grid_view()->template begin<0>(); entity_it != entity_end_it; ++entity_it) {
      const auto& entity = *entity_it;
      // * as the derived type
      BaseFunctionSetType d_base_function_set = space_->base_function_set(entity);
      const D_BackendType& d_backend          = d_base_function_set.backend();
      size_t d_order = d_base_function_set.order();
      if (d_order != SpaceType::polOrder)
        DUNE_THROW_COLORFULLY(Exceptions::internal_error, d_order << " vs. " << SpaceType::polOrder);
      //   the size has already been checked in fulfills_interface() above
      // * as the interface
      InterfaceType& i_base_function_set = static_cast<InterfaceType&>(d_base_function_set);
      const I_BackendType& i_backend = i_base_function_set.backend();
      if (&d_backend != &i_backend)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
      size_t i_order = i_base_function_set.order();
      if (i_order != d_order)
        DUNE_THROW_COLORFULLY(Exceptions::CRTP_check_failed, "");
    } // walk the grid
  } // ... basefunctionset_fulfills_interface()

protected:
  ProviderType grid_provider_;
  std::unique_ptr<const SpaceType> space_;
}; // struct SpaceTestBase
