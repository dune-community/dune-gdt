#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH

#include <dune/gdt/discretefunction/default.hh>

#include "l2-local.hh"
#include "l2-global.hh"

namespace Dune {
namespace GDT {

// forward
template< class GridViewImp, class FieldImp = double >
class L2ProjectionOperator;


namespace internal {


template< class GridViewType, class SourceType, class RangeType >
class L2ProjectionLocalizableOperatorTraits
{
  static_assert(is_discrete_function< RangeType >::value, "");

  template< class G, class S, class R, bool c = true >
  struct Helper {
    typedef L2GlobalProjectionLocalizableOperator< G, S, R > type;
  };

  template< class G, class S, class R >
  struct Helper< G, S, R, false > {
    typedef L2LocalProjectionLocalizableOperator< G, S, R > type;
  };

public:
  typedef typename Helper< GridViewType, SourceType, RangeType, RangeType::SpaceType::continuous >::type BaseType;
}; // class L2ProjectionLocalizableOperatorTraits


} // namespace internal


template< class GridViewImp,
          class SourceImp,
          class RangeImp >
class L2ProjectionLocalizableOperator
  : public internal::L2ProjectionLocalizableOperatorTraits< GridViewImp, SourceImp, RangeImp >::BaseType
{
  typedef typename internal::L2ProjectionLocalizableOperatorTraits< GridViewImp, SourceImp, RangeImp >::BaseType BaseType;
public:
  template< class ...Args >
  explicit L2ProjectionLocalizableOperator(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  L2ProjectionLocalizableOperator< GridViewType, SourceType,
                                                                   DiscreteFunction< SpaceType, VectorType > >
                                            > >::type
make_l2_projection_localizable_operator(const GridViewType& grid_view,
                                        const SourceType& source,
                                        DiscreteFunction< SpaceType, VectorType >& range,
                                        const size_t over_integrate = 0)
{
  return DSC::make_unique<
      L2ProjectionLocalizableOperator< GridViewType, SourceType,
                                       DiscreteFunction< SpaceType, VectorType > > >(over_integrate,
                                                                                     grid_view,
                                                                                     source,
                                                                                     range);
} // ... make_l2_projection_localizable_operator(...)

template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  L2ProjectionLocalizableOperator< typename SpaceType::GridViewType,
                                                                   SourceType,
                                                                   DiscreteFunction< SpaceType, VectorType > >
                                            > >::type
make_l2_projection_localizable_operator(const SourceType& source,
                                        DiscreteFunction< SpaceType, VectorType >& range,
                                        const size_t over_integrate = 0)
{
  return DSC::make_unique<
      L2ProjectionLocalizableOperator< typename SpaceType::GridViewType,
                                       SourceType,
                                       DiscreteFunction< SpaceType, VectorType > > >(over_integrate,
                                                                                     range.space().grid_view(),
                                                                                     source,
                                                                                     range);
} // ... make_l2_projection_localizable_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
