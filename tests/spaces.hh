// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <type_traits>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>

#include <dune/gdt/space/tools.hh>
#include <dune/gdt/space/interface.hh>
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
typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, false>::Type AluConform2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, true>::Type AluConform2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluConform2dGridType, true>::Type AluConform2dLevelGridViewType;
typedef Dune::ALUSimplexGrid<2, 2> AluSimplex2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, false>::Type AluSimplex2dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex2dGridType, true>::Type AluSimplex2dLevelGridViewType;
typedef Dune::ALUSimplexGrid<3, 3> AluSimplex3dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLeafGridPartType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, false>::Type AluSimplex3dLevelGridPartType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLeafGridViewType;
typedef
    typename Dune::GDT::SpaceTools::LevelGridPartView<AluSimplex3dGridType, true>::Type AluSimplex3dLevelGridViewType;
typedef Dune::ALUCubeGrid<3, 3> AluCube3dGridType;
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

public:
  ~SpaceTestBase()
  {
  }

  SpaceTestBase()
  {
    using namespace Dune;
    using namespace GDT;
    Stuff::GridProviderCube<GridType> grid_provider(0.0, 1.0, 4u);
    grid_                     = grid_provider.grid();
    const auto grid_part_view = SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_);
    space_                    = std::unique_ptr<SpaceType>(new SpaceType(grid_part_view));
  }

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
    // walk the grid
    const auto d_entity_it_end = d_grid_view->template end<0>();
    for (auto d_entity_it = d_grid_view->template begin<0>(); d_entity_it != d_entity_it_end; ++d_entity_it) {
      const D_EntityType& d_entity = *d_entity_it;
      D_BaseFunctionSetType DUNE_UNUSED(d_base_function_set) = space_->base_function_set(d_entity);
    }
    D_PatternType d_pattern                        = space_->compute_pattern();
    D_PatternType d_pattern_view                   = space_->compute_pattern(*d_grid_view);
    D_PatternType d_pattern_other                  = space_->compute_pattern(*space_);
    D_PatternType d_pattern_view_other             = space_->compute_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_volume                 = space_->compute_volume_pattern();
    D_PatternType d_pattern_volume_view            = space_->compute_volume_pattern(*d_grid_view);
    D_PatternType d_pattern_volume_other           = space_->compute_volume_pattern(*space_);
    D_PatternType d_pattern_volume_view_other      = space_->compute_volume_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_face_volume            = space_->compute_face_and_volume_pattern();
    D_PatternType d_pattern_face_volume_view       = space_->compute_face_and_volume_pattern(*d_grid_view);
    D_PatternType d_pattern_face_volume_other      = space_->compute_face_and_volume_pattern(*space_);
    D_PatternType d_pattern_face_volume_view_other = space_->compute_face_and_volume_pattern(*d_grid_view, *space_);
    D_PatternType d_pattern_face                   = space_->compute_face_pattern();
    D_PatternType d_pattern_face_view              = space_->compute_face_pattern(*d_grid_view);
    D_PatternType d_pattern_face_other             = space_->compute_face_pattern(*space_);
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
    // walk the grid
    const auto i_entity_it_end = i_grid_view->template end<0>();
    for (auto i_entity_it = d_grid_view->template begin<0>(); i_entity_it != i_entity_it_end; ++i_entity_it) {
      const I_EntityType& i_entity = *i_entity_it;
      I_BaseFunctionSetType DUNE_UNUSED(i_base_function_set) = i_space.base_function_set(i_entity);
    }
    I_PatternType i_pattern                        = i_space.compute_pattern();
    I_PatternType i_pattern_view                   = i_space.compute_pattern(*i_grid_view);
    I_PatternType i_pattern_other                  = i_space.compute_pattern(i_space);
    I_PatternType i_pattern_view_other             = i_space.compute_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_volume                 = i_space.compute_volume_pattern();
    I_PatternType i_pattern_volume_view            = i_space.compute_volume_pattern(*i_grid_view);
    I_PatternType i_pattern_volume_other           = i_space.compute_volume_pattern(i_space);
    I_PatternType i_pattern_volume_view_other      = i_space.compute_volume_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_face_volume            = i_space.compute_face_and_volume_pattern();
    I_PatternType i_pattern_face_volume_view       = i_space.compute_face_and_volume_pattern(*i_grid_view);
    I_PatternType i_pattern_face_volume_other      = i_space.compute_face_and_volume_pattern(i_space);
    I_PatternType i_pattern_face_volume_view_other = i_space.compute_face_and_volume_pattern(*i_grid_view, i_space);
    I_PatternType i_pattern_face                   = i_space.compute_face_pattern();
    I_PatternType i_pattern_face_view              = i_space.compute_face_pattern(*i_grid_view);
    I_PatternType i_pattern_face_other             = i_space.compute_face_pattern(i_space);
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
  } // ... fulfills_interface()

protected:
  std::shared_ptr<GridType> grid_;
  std::unique_ptr<const SpaceType> space_;
}; // struct SpaceTestBase
