// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH
#define DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>

#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/projections/l2-local.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class L2LocalProlongationOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class L2LocalProlongationOperatorTraits
{
public:
  typedef L2LocalProlongationOperator<GridLayerImp, FieldImp> derived_type;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \brief Carries out a prolongation (in a localized manner) using a local L2 projection.
 *
 *        This is done by reinterpreting the source on the range grid layer and applying a
 *        L2LocalProjectionLocalizableOperator.
 */
template <class GridLayerImp, class SourceImp, class RangeImp>
class L2LocalProlongationLocalizableOperator
  : XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>>
  , public L2LocalProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
{
  static_assert(is_const_discrete_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>> SourceStorage;
  typedef L2LocalProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
      BaseOperatorType;

public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridLayerType;
  using typename BaseOperatorType::RangeType;

  L2LocalProlongationLocalizableOperator(const size_t over_integrate,
                                         GridLayerType grd_vw,
                                         const SourceType& src,
                                         RangeType& rng,
                                         const XT::Common::Parameter& param = {})
    : SourceStorage(new ReinterpretDiscreteFunction<SourceImp>(src))
    , BaseOperatorType(over_integrate, grd_vw, SourceStorage::access(), rng, param)
  {}

  L2LocalProlongationLocalizableOperator(GridLayerType grd_vw,
                                         const SourceType& src,
                                         RangeType& rng,
                                         const XT::Common::Parameter& param = {})
    : L2LocalProlongationLocalizableOperator(0, grd_vw, src, rng, param)
  {}

  ///! Calls L2LocalProjectionLocalizableOperator::apply and gives a meaningful error message.
  void apply()
  {
    try {
      BaseOperatorType::apply();
    } catch (XT::Functions::Exceptions::reinterpretation_error& ee) {
      DUNE_THROW(prolongation_error,
                 "This prolongation (using a local L2 projection) failed, because the source could not be reinterpreted"
                     << " on the given grid layer!\n"
                     << "This was the original error:\n\n"
                     << ee.what());
    }
  } // ... apply(...)
}; // class L2LocalProlongationLocalizableOperator


template <class GridLayerType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<
    XT::Grid::is_layer<GridLayerType>::value && is_space<SourceSpaceType>::value
        && XT::LA::is_vector<SourceVectorType>::value && is_space<RangeSpaceType>::value
        && XT::LA::is_vector<RangeVectorType>::value,
    std::unique_ptr<L2LocalProlongationLocalizableOperator<GridLayerType,
                                                           ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                           DiscreteFunction<RangeSpaceType, RangeVectorType>>>>::type
make_local_l2_prolongation_localizable_operator(const GridLayerType& grid_layer,
                                                const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
                                                DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
                                                const size_t over_integrate = 0,
                                                const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<
      L2LocalProlongationLocalizableOperator<GridLayerType,
                                             ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                             DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
      over_integrate, grid_layer, source, range, param);
} // ... make_local_l2_prolongation_localizable_operator(...)

template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
typename std::enable_if<
    is_space<SourceSpaceType>::value && XT::LA::is_vector<SourceVectorType>::value && is_space<RangeSpaceType>::value
        && XT::LA::is_vector<RangeVectorType>::value,
    std::unique_ptr<L2LocalProlongationLocalizableOperator<typename RangeSpaceType::GridLayerType,
                                                           ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                           DiscreteFunction<RangeSpaceType, RangeVectorType>>>>::type
make_local_l2_prolongation_localizable_operator(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
                                                DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
                                                const size_t over_integrate = 0,
                                                const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<
      L2LocalProlongationLocalizableOperator<typename RangeSpaceType::GridLayerType,
                                             ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                             DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
      over_integrate, range.space().grid_layer(), source, range, param);
} // ... make_local_l2_prolongation_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2LocalProlongationOperator
  : public OperatorInterface<internal::L2LocalProlongationOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2LocalProlongationOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2LocalProlongationOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2LocalProlongationOperator(const size_t over_integrate, GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {}

  L2LocalProlongationOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
  {}

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source,
             DiscreteFunction<RS, RV>& range,
             const XT::Common::Parameter& param = {}) const
  {
    L2LocalProlongationLocalizableOperator<GridLayerType, ConstDiscreteFunction<SS, SV>, DiscreteFunction<RS, RV>> op(
        over_integrate_, grid_layer_, source, range, param);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType
  apply2(const RangeType& /*range*/, const SourceType& /*source*/, const XT::Common::Parameter& /*param*/ = {}) const
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
  const size_t over_integrate_;
}; // class L2LocalProlongationOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2LocalProlongationOperator<GridLayerType>>>::type
make_local_l2_prolongation_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2LocalProlongationOperator<GridLayerType>>(over_integrate, grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_LOCAL_HH
