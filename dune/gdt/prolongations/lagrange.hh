// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
#define DUNE_GDT_PROLONGATIONS_LAGRANGE_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections/lagrange.hh>


namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class LagrangeProlongationOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class LagrangeProlongationOperatorTraits
{
public:
  typedef LagrangeProlongationOperator<GridLayerImp, FieldImp> derived_type;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \brief Carries out a prolongation (in a localized manner) using a lagrange projection.
 *
 *        This is done by reinterpreting the source on the range grid layer and applying a
 *        LagrangeProjectionLocalizableOperator.
 */
template <class GridLayerImp, class SourceImp, class RangeImp>
class LagrangeProlongationLocalizableOperator
    : XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>>,
      public LagrangeProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
{
  static_assert(is_const_discrete_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>> SourceStorage;
  typedef LagrangeProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
      BaseOperatorType;

public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridLayerType;
  using typename BaseOperatorType::RangeType;

  LagrangeProlongationLocalizableOperator(GridLayerType grd_vw, const SourceType& src, RangeType& rng, const XT::Common::Parameter& param = {})
    : SourceStorage(new ReinterpretDiscreteFunction<SourceImp>(src))
    , BaseOperatorType(grd_vw, SourceStorage::access(), rng, param)
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
                     << " on the given grid layer!\n"
                     << "This was the original error:\n\n"
                     << ee.what());
    }
  } // ... apply(...)
}; // class LagrangeProlongationLocalizableOperator


template <class GridLayerType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<LagrangeProlongationLocalizableOperator<GridLayerType,
                                                                                ConstDiscreteFunction<SourceSpaceType,
                                                                                                      SourceVectorType>,
                                                                                DiscreteFunction<RangeSpaceType,
                                                                                                 RangeVectorType>>>>::
    type
    make_lagrange_prolongation_localizable_operator(
        const GridLayerType& grid_layer,
        const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
    DiscreteFunction<RangeSpaceType, RangeVectorType>& range, const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::
      make_unique<LagrangeProlongationLocalizableOperator<GridLayerType,
                                                          ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                          DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
          grid_layer, source, range, param);
} // ... make_lagrange_prolongation_localizable_operator(...)

template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
std::unique_ptr<LagrangeProlongationLocalizableOperator<typename RangeSpaceType::GridLayerType,
                                                        ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                        DiscreteFunction<RangeSpaceType, RangeVectorType>>>
make_lagrange_prolongation_localizable_operator(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
                                                DiscreteFunction<RangeSpaceType, RangeVectorType>& rangeDiscreteFunction<RangeSpaceType, RangeVectorType>& range, const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<LagrangeProlongationLocalizableOperator<
      typename RangeSpaceType::GridLayerType,
      ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
      DiscreteFunction<RangeSpaceType, RangeVectorType>>>(range.space().grid_layer(), source, range, param);
} // ... make_lagrange_prolongation_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class LagrangeProlongationOperator
    : public OperatorInterface<internal::LagrangeProlongationOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::LagrangeProlongationOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::LagrangeProlongationOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  LagrangeProlongationOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
  {
  }

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range, const XT::Common::Parameter& param = {}) const
  {
    LagrangeProlongationLocalizableOperator<GridLayerType, ConstDiscreteFunction<SS, SV>, DiscreteFunction<RS, RV>> op(
        grid_layer_, source, range, param);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/, const XT::Common::Parameter& /*param*/ = {}) const
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
  GridLayerType grid_layer_;
}; // class LagrangeProlongationOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<LagrangeProlongationOperator<GridLayerType>>>::type
make_lagrange_prolongation_operator(const GridLayerType& grid_layer)
{
  return Dune::XT::Common::make_unique<LagrangeProlongationOperator<GridLayerType>>(grid_layer);
}


template <class GridLayerType, class SS, class SV, class RS, class RV>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value, void>::type prolong_lagrange(
    const GridLayerType& grid_layer, const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range, const XT::Common::Parameter& param = {})
{
  make_lagrange_prolongation_operator(grid_layer)->apply(source, range, param);
}

template <class SS, class SV, class RS, class RV>
void prolong_lagrange(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range, const XT::Common::Parameter& param = {})
{
  make_lagrange_prolongation_operator(range.space().grid_layer())->apply(source, range, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_LAGRANGE_HH
