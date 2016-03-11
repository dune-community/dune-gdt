// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH
#define DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH

#include <dune/stuff/common/memory.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/projections/l2-local.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Carries out a prolongation (in a localized manner) using a local L2 projection.
 *
 *        This is done by reinterpreting the source on the range grid view and applying a
 *        L2LocalProjectionLocalizableOperator.
 */
template< class GridViewImp, class SourceImp, class RangeImp >
class L2LocalProlongationLocalizableOperator
  : Stuff::Common::ConstStorageProvider< ReinterpretDiscreteFunction< SourceImp > >
  , public L2LocalProjectionLocalizableOperator< GridViewImp, ReinterpretDiscreteFunction< SourceImp >, RangeImp >
{
  typedef Stuff::Common::ConstStorageProvider< ReinterpretDiscreteFunction< SourceImp > > SourceStorage;
  typedef L2LocalProjectionLocalizableOperator
      < GridViewImp, ReinterpretDiscreteFunction< SourceImp >, RangeImp >                 BaseOperatorType;
public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridViewType;
  using typename BaseOperatorType::RangeType;

  L2LocalProlongationLocalizableOperator(const size_t over_integrate,
                                         GridViewType grd_vw,
                                         const SourceType& src,
                                         RangeType& rng)
    : SourceStorage(new ReinterpretDiscreteFunction< SourceType >(src))
    , BaseOperatorType(over_integrate, grd_vw, SourceStorage::access(), rng)
  {}

  L2LocalProlongationLocalizableOperator(GridViewType grd_vw, const SourceType& src, RangeType& rng)
    : L2LocalProlongationLocalizableOperator(0, grd_vw, src, rng)
  {}

  ///! Calls L2LocalProjectionLocalizableOperator::apply and gives a meaningful error message.
  void apply()
  {
    try {
      BaseOperatorType::apply();
    } catch (Stuff::Exceptions::reinterpretation_error& ee) {
      DUNE_THROW(Exceptions::prolongation_error,
                 "This prolongation (using a local L2 projection) failed, because the source could not be reinterpreted"
                 << " on the given grid view!\n"
                 << "This was the original error:\n\n" << ee.what());
    }
  } // ... apply(...)
}; // class L2LocalProlongationLocalizableOperator


template< class GridViewType, class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  L2LocalProlongationLocalizableOperator< GridViewType, SourceType,
                                                                          DiscreteFunction< SpaceType, VectorType > >
                                            > >::type
make_local_l2_prolongation_localizable_operator(const GridViewType& grid_view,
                                                const SourceType& source,
                                                DiscreteFunction< SpaceType, VectorType >& range,
                                                const size_t over_integrate = 0)
{
  return DSC::make_unique<
      L2LocalProlongationLocalizableOperator< GridViewType, SourceType,
                                              DiscreteFunction< SpaceType, VectorType > > >(over_integrate,
                                                                                            grid_view,
                                                                                            source,
                                                                                            range);
} // ... make_local_l2_prolongation_localizable_operator(...)

template< class SourceType, class SpaceType, class VectorType >
    typename std::enable_if<    Stuff::is_localizable_function< SourceType >::value
                             && is_space< SpaceType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                           , std::unique_ptr<
                                  L2LocalProlongationLocalizableOperator< typename SpaceType::GridViewType,
                                                                          SourceType,
                                                                          DiscreteFunction< SpaceType, VectorType > >
                                            > >::type
make_local_l2_prolongation_localizable_operator(const SourceType& source,
                                                DiscreteFunction< SpaceType, VectorType >& range,
                                                const size_t over_integrate = 0)
{
  return DSC::make_unique<
      L2LocalProlongationLocalizableOperator< typename SpaceType::GridViewType,
                                            SourceType,
                                            DiscreteFunction< SpaceType, VectorType > > >(over_integrate,
                                                                                          range.space().grid_view(),
                                                                                          source,
                                                                                          range);
} // ... make_local_l2_prolongation_localizable_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH
