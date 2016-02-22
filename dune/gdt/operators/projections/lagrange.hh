// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_LAGRANGE_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_LAGRANGE_HH

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/lagrange-projection.hh>
#include <dune/gdt/spaces/interface.hh>

#include "../default.hh"
#include "../interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \todo Add a check if the range space provides Lagrange points (after implementing the appropriate interface/mixin
 *       for those spaces).
 * \note Do we have to set all range DoFs to infinity here?
 */
template< class GridViewImp, class SourceImp, class RangeImp >
class LagrangeProjectionLocalizableOperator
  : public LocalizableOperatorDefault< GridViewImp, SourceImp, RangeImp >
{
  typedef LocalizableOperatorDefault< GridViewImp, SourceImp, RangeImp > BaseType;
public:
  template< class ...Args >
  explicit LagrangeProjectionLocalizableOperator(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_operator_()
  {
    this->add(local_operator_);
  }

private:
  const LocalLagrangeProjectionOperator local_operator_;
}; // class LagrangeProjectionLocalizableOperator


template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  LagrangeProjectionLocalizableOperator< GridViewType, SourceType,
                                                                         DiscreteFunction< SpaceType, VectorType > > >
                           >::type
make_lagrange_projection_localizable_operator(const GridViewType& grid_view,
                                              const SourceType& source,
                                              DiscreteFunction< SpaceType, VectorType >& range)
{
  return DSC::make_unique< LagrangeProjectionLocalizableOperator
      < GridViewType, SourceType, DiscreteFunction< SpaceType, VectorType > > >(
            grid_view, source, range);
}


template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  LagrangeProjectionLocalizableOperator< typename SpaceType::GridViewType, SourceType,
                                                                         DiscreteFunction< SpaceType, VectorType > > >
                           >::type
make_lagrange_projection_localizable_operator(const SourceType& source,
                                              DiscreteFunction< SpaceType, VectorType >& range)
{
  return DSC::make_unique< LagrangeProjectionLocalizableOperator
      < typename SpaceType::GridViewType, SourceType, DiscreteFunction< SpaceType, VectorType > > >(
          range.space().grid_view(), source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_LAGRANGE_HH
