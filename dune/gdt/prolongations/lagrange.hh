// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
#define DUNE_GDT_PROLONGATIONS_LAGRANGE_HH

#include <dune/xt/common/memory.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections/lagrange.hh>


namespace Dune {
namespace GDT {


// forward
template <class GridViewImp, class FieldImp = double>
class LagrangeProlongationOperator;


namespace internal {


template <class GridViewImp, class FieldImp>
class LagrangeProlongationOperatorTraits
{
public:
  typedef LagrangeProlongationOperator<GridViewImp, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \brief Carries out a prolongation (in a localized manner) using a lagrange projection.
 *
 *        This is done by reinterpreting the source on the range grid view and applying a
 *        LagrangeProjectionLocalizableOperator.
 */
template <class GridViewImp, class SourceImp, class RangeImp>
class LagrangeProlongationLocalizableOperator
    : XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>>,
      public LagrangeProjectionLocalizableOperator<GridViewImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
{
  static_assert(is_const_discrete_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>> SourceStorage;
  typedef LagrangeProjectionLocalizableOperator<GridViewImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
      BaseOperatorType;

public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridViewType;
  using typename BaseOperatorType::RangeType;

  LagrangeProlongationLocalizableOperator(GridViewType grd_vw, const SourceType& src, RangeType& rng)
    : SourceStorage(new ReinterpretDiscreteFunction<SourceImp>(src))
    , BaseOperatorType(grd_vw, SourceStorage::access(), rng)
  {
  }

  ///! Calls LagrangeProjectionLocalizableOperator::apply and gives a meaningful error message.
  void apply()
  {
    try {
      BaseOperatorType::apply();
    } catch (XT::Common::Exceptions::reinterpretation_error& ee) {
      DUNE_THROW(prolongation_error,
                 "This prolongation (using a lagrange projection) failed, because the source could not be reinterpreted"
                     << " on the given grid view!\n"
                     << "This was the original error:\n\n"
                     << ee.what());
    }
  } // ... apply(...)
}; // class LagrangeProlongationLocalizableOperator


template <class GridViewType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<LagrangeProlongationLocalizableOperator<GridViewType,
                                                                                ConstDiscreteFunction<SourceSpaceType,
                                                                                                      SourceVectorType>,
                                                                                DiscreteFunction<RangeSpaceType,
                                                                                                 RangeVectorType>>>>::
    type
    make_lagrange_prolongation_localizable_operator(
        const GridViewType& grid_view,
        const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range)
{
  return Dune::XT::Common::
      make_unique<LagrangeProlongationLocalizableOperator<GridViewType,
                                                          ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                          DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
          grid_view, source, range);
} // ... make_lagrange_prolongation_localizable_operator(...)

template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
std::unique_ptr<LagrangeProlongationLocalizableOperator<typename RangeSpaceType::GridViewType,
                                                        ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                        DiscreteFunction<RangeSpaceType, RangeVectorType>>>
make_lagrange_prolongation_localizable_operator(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
                                                DiscreteFunction<RangeSpaceType, RangeVectorType>& range)
{
  return Dune::XT::Common::make_unique<LagrangeProlongationLocalizableOperator<
      typename RangeSpaceType::GridViewType,
      ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
      DiscreteFunction<RangeSpaceType, RangeVectorType>>>(range.space().grid_view(), source, range);
} // ... make_lagrange_prolongation_localizable_operator(...)


template <class GridViewImp, class FieldImp>
class LagrangeProlongationOperator
    : public OperatorInterface<internal::LagrangeProlongationOperatorTraits<GridViewImp, FieldImp>>
{
  typedef OperatorInterface<internal::LagrangeProlongationOperatorTraits<GridViewImp, FieldImp>> BaseType;

public:
  typedef internal::LagrangeProlongationOperatorTraits<GridViewImp, FieldImp> Traits;
  typedef GridViewImp GridViewType;
  using typename BaseType::FieldType;

private:
  typedef typename XT::Grid::Entity<GridViewType>::Type E;
  typedef typename GridViewType::ctype D;
  static const size_t d = GridViewType::dimension;

public:
  LagrangeProlongationOperator(GridViewType grid_view)
    : grid_view_(grid_view)
  {
  }

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range) const
  {
    LagrangeProlongationLocalizableOperator<GridViewType, ConstDiscreteFunction<SS, SV>, DiscreteFunction<RS, RV>> op(
        grid_view_, source, range);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template <class RangeType, class SourceType>
  void
  apply_inverse(const RangeType& /*range*/, SourceType& /*source*/, const XT::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

private:
  GridViewType grid_view_;
}; // class LagrangeProlongationOperator


template <class GridViewType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<LagrangeProlongationOperator<GridViewType>>>::type
make_lagrange_prolongation_operator(const GridViewType& grid_view)
{
  return Dune::XT::Common::make_unique<LagrangeProlongationOperator<GridViewType>>(grid_view);
}


template <class GridViewType, class SS, class SV, class RS, class RV>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value, void>::type prolong_lagrange(
    const GridViewType& grid_view, const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range)
{
  make_lagrange_prolongation_operator(grid_view)->apply(source, range);
}

template <class SS, class SV, class RS, class RV>
void prolong_lagrange(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range)
{
  make_lagrange_prolongation_operator(range.space().grid_view())->apply(source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
