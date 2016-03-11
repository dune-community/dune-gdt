// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
#define DUNE_GDT_PROLONGATIONS_LAGRANGE_HH

#include <dune/stuff/common/memory.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/projections/lagrange.hh>


namespace Dune {
namespace GDT {


/**
 * \brief Carries out a prolongation (in a localized manner) using a lagrange projection.
 *
 *        This is done by reinterpreting the source on the range grid view and applying a
 *        LagrangeProjectionLocalizableOperator.
 */
template< class GridViewImp, class SourceImp, class RangeImp >
class LagrangeProlongationLocalizableOperator
  : Stuff::Common::ConstStorageProvider< ReinterpretDiscreteFunction< SourceImp > >
  , public LagrangeProjectionLocalizableOperator< GridViewImp, ReinterpretDiscreteFunction< SourceImp >, RangeImp >
{
  static_assert(is_const_discrete_function< SourceImp >::value, "");
  static_assert(is_discrete_function< RangeImp >::value, "");
  typedef Stuff::Common::ConstStorageProvider< ReinterpretDiscreteFunction< SourceImp > > SourceStorage;
  typedef LagrangeProjectionLocalizableOperator
      < GridViewImp, ReinterpretDiscreteFunction< SourceImp >, RangeImp >                 BaseOperatorType;
public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridViewType;
  using typename BaseOperatorType::RangeType;

  LagrangeProlongationLocalizableOperator(GridViewType grd_vw, const SourceType& src, RangeType& rng)
    : SourceStorage(new ReinterpretDiscreteFunction< SourceImp >(src))
    , BaseOperatorType(grd_vw, SourceStorage::access(), rng)
  {}

  ///! Calls LagrangeProjectionLocalizableOperator::apply and gives a meaningful error message.
  void apply()
  {
    try {
      BaseOperatorType::apply();
    } catch (Stuff::Exceptions::reinterpretation_error& ee) {
      DUNE_THROW(Exceptions::prolongation_error,
                 "This prolongation (using a lagrange projection) failed, because the source could not be reinterpreted"
                 << " on the given grid view!\n"
                 << "This was the original error:\n\n" << ee.what());
    }
  } // ... apply(...)
}; // class LagrangeProlongationLocalizableOperator


template< class GridViewType, class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType >
    typename std::enable_if< Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr<
                                  LagrangeProlongationLocalizableOperator< GridViewType,
                                                                           ConstDiscreteFunction< SourceSpaceType, SourceVectorType >,
                                                                           DiscreteFunction< RangeSpaceType, RangeVectorType > >
                                            > >::type
make_lagrange_prolongation_localizable_operator(const GridViewType& grid_view,
                                                const ConstDiscreteFunction< SourceSpaceType, SourceVectorType >& source,
                                                DiscreteFunction< RangeSpaceType, RangeVectorType >& range)
{
  return DSC::make_unique<
      LagrangeProlongationLocalizableOperator< GridViewType,
                                               ConstDiscreteFunction< SourceSpaceType, SourceVectorType >,
                                               DiscreteFunction< RangeSpaceType, RangeVectorType > > >(grid_view,
                                                                                                       source,
                                                                                                       range);
} // ... make_lagrange_prolongation_localizable_operator(...)

template< class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType >
    std::unique_ptr< LagrangeProlongationLocalizableOperator< typename RangeSpaceType::GridViewType,
                                                              ConstDiscreteFunction< SourceSpaceType, SourceVectorType >,
                                                              DiscreteFunction< RangeSpaceType, RangeVectorType > > >
make_lagrange_prolongation_localizable_operator(const ConstDiscreteFunction< SourceSpaceType, SourceVectorType >& source,
                                                DiscreteFunction< RangeSpaceType, RangeVectorType >& range)
{
  return DSC::make_unique<
      LagrangeProlongationLocalizableOperator< typename RangeSpaceType::GridViewType,
                                               ConstDiscreteFunction< SourceSpaceType, SourceVectorType >,
                                               DiscreteFunction< RangeSpaceType, RangeVectorType > > >(range.space().grid_view(),
                                                                                                       source,
                                                                                                       range);
} // ... make_lagrange_prolongation_localizable_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
