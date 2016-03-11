// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PROJECTIONS_LAGRANGE_HH
#define DUNE_GDT_PROJECTIONS_LAGRANGE_HH

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/lagrange-projection.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/operators/interfaces.hh>

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


// forward
template< class GridViewImp, class FieldImp = double >
class LagrangeProjectionOperator;


namespace internal {


template< class GridViewImp, class FieldImp >
class LagrangeProjectionOperatorTraits
{
public:
  typedef LagrangeProjectionOperator< GridViewImp, FieldImp > derived_type;
  typedef FieldImp                                           FieldType;
};


} // namespace internal


template< class GridViewImp, class FieldImp >
class LagrangeProjectionOperator
  : public OperatorInterface< internal::LagrangeProjectionOperatorTraits< GridViewImp, FieldImp > >
{
  typedef OperatorInterface< internal::LagrangeProjectionOperatorTraits< GridViewImp, FieldImp > > BaseType;
public:
  typedef internal::LagrangeProjectionOperatorTraits< GridViewImp, FieldImp > Traits;
  typedef GridViewImp GridViewType;
  using typename BaseType::FieldType;
private:
  typedef typename Stuff::Grid::Entity< GridViewType >::Type E;
  typedef typename GridViewType::ctype D;
  static const size_t d = GridViewType::dimension;

public:
  LagrangeProjectionOperator(GridViewType grid_view)
    : grid_view_(grid_view)
  {}

  template< class R, size_t r, size_t rC, class S, class V >
  void apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
             DiscreteFunction< S, V >& range) const
  {
    make_lagrange_projection_localizable_operator(grid_view_, source, range)->apply();
  }

  template< class RangeType, class SourceType >
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template< class RangeType, class SourceType >
  void apply_inverse(const RangeType& /*range*/,
                     SourceType& /*source*/,
                     const Stuff::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  std::vector< std::string > invert_options() const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  Stuff::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

private:
  GridViewType grid_view_;
}; // class LagrangeProjectionOperator


template< class GridViewType >
    typename std::enable_if< Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< LagrangeProjectionOperator< GridViewType > >
                           >::type
make_lagrange_projection_operator(const GridViewType& grid_view)
{
  return DSC::make_unique< LagrangeProjectionOperator< GridViewType > >(grid_view);
}


template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project_lagrange(const GridViewType& grid_view,
                 const SourceType& source,
                 DiscreteFunction< SpaceType, VectorType >& range)
{
  make_lagrange_projection_operator(grid_view)->apply(source, range);
}


template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , void >::type
project_lagrange(const SourceType& source, DiscreteFunction< SpaceType, VectorType >& range)
{
  make_lagrange_projection_operator(range.space().grid_view())->apply(source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_LAGRANGE_HH
