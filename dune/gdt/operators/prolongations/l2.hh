// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROLONGATIONS_L2_HH
#define DUNE_GDT_OPERATORS_PROLONGATIONS_L2_HH

#include <dune/stuff/common/memory.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/operators/projections/l2-local.hh>

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
                                         GridViewType grid_view,
                                         const SourceType& source,
                                         RangeType& range)
    : SourceStorage(new ReinterpretDiscreteFunction< SourceType >(source))
    , BaseOperatorType(over_integrate, grid_view, SourceStorage::access(), range)
  {}

  L2LocalProlongationLocalizableOperator(GridViewType grid_view, const SourceType& source, RangeType& range)
    : SourceStorage(new ReinterpretDiscreteFunction< SourceType >(source))
    , BaseOperatorType(grid_view, SourceStorage::access(), range)
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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROLONGATIONS_L2_HH
