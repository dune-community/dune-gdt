// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROLONGATIONS_L2_GLOBAL_HH
#define DUNE_GDT_PROLONGATIONS_L2_GLOBAL_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/projections/l2-global.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class L2GlobalProlongationOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class L2GlobalProlongationOperatorTraits
{
public:
  typedef L2GlobalProlongationOperator<GridLayerImp, FieldImp> derived_type;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class SourceImp, class RangeImp>
class L2GlobalProlongationLocalizableOperator
    : XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>>,
      public L2GlobalProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
{
  static_assert(is_const_discrete_function<SourceImp>::value, "");
  static_assert(is_discrete_function<RangeImp>::value, "");
  typedef XT::Common::ConstStorageProvider<ReinterpretDiscreteFunction<SourceImp>> SourceStorage;
  typedef L2GlobalProjectionLocalizableOperator<GridLayerImp, ReinterpretDiscreteFunction<SourceImp>, RangeImp>
      BaseOperatorType;

public:
  typedef SourceImp SourceType;
  using typename BaseOperatorType::GridLayerType;
  using typename BaseOperatorType::RangeType;

  L2GlobalProlongationLocalizableOperator(const size_t over_integrate,
                                          GridLayerType grd_vw,
                                          const SourceType& src,
                                          RangeType& rng,
                                          const XT::Common::Parameter& param = {})
    : SourceStorage(new ReinterpretDiscreteFunction<SourceImp>(src))
    , BaseOperatorType(over_integrate, grd_vw, SourceStorage::access(), rng, param)
  {
  }

  L2GlobalProlongationLocalizableOperator(GridLayerType grd_vw,
                                          const SourceType& src,
                                          RangeType& rng,
                                          const XT::Common::Parameter& param = {})
    : L2GlobalProlongationLocalizableOperator(0, grd_vw, src, rng, param)
  {
  }

  ///! Calls L2GlobalProjectionLocalizableOperator::apply and gives a meaningful error message.
  void apply()
  {
    try {
      BaseOperatorType::apply();
    } catch (XT::Common::Exceptions::reinterpretation_error& ee) {
      DUNE_THROW(prolongation_error,
                 "This prolongation (using a global L2 projection) failed, because the source could not be "
                     << "reinterpreted on the given grid layer!\n"
                     << "This was the original error:\n\n"
                     << ee.what());
    }
  } // ... apply(...)
}; // class L2GlobalProlongationLocalizableOperator


template <class GridLayerType,
          class SourceSpaceType,
          class SourceVectorType,
          class RangeSpaceType,
          class RangeVectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value && is_space<SourceSpaceType>::value
                            && XT::LA::is_vector<SourceVectorType>::value
                            && is_space<RangeSpaceType>::value
                            && XT::LA::is_vector<RangeVectorType>::value,
                        std::unique_ptr<L2GlobalProlongationLocalizableOperator<GridLayerType,
                                                                                ConstDiscreteFunction<SourceSpaceType,
                                                                                                      SourceVectorType>,
                                                                                DiscreteFunction<RangeSpaceType,
                                                                                                 RangeVectorType>>>>::
    type
    make_global_l2_prolongation_localizable_operator(
        const GridLayerType& grid_layer,
        const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
        DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
        const size_t over_integrate = 0,
        const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::
      make_unique<L2GlobalProlongationLocalizableOperator<GridLayerType,
                                                          ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                          DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
          over_integrate, grid_layer, source, range, param);
} // ... make_global_l2_prolongation_localizable_operator(...)

template <class SourceSpaceType, class SourceVectorType, class RangeSpaceType, class RangeVectorType>
typename std::enable_if<is_space<SourceSpaceType>::value && XT::LA::is_vector<SourceVectorType>::value
                            && is_space<RangeSpaceType>::value
                            && XT::LA::is_vector<RangeVectorType>::value,
                        std::unique_ptr<L2GlobalProlongationLocalizableOperator<
                            typename RangeSpaceType::GridLayerType,
                            ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                            DiscreteFunction<RangeSpaceType, RangeVectorType>>>>::type
make_global_l2_prolongation_localizable_operator(const ConstDiscreteFunction<SourceSpaceType, SourceVectorType>& source,
                                                 DiscreteFunction<RangeSpaceType, RangeVectorType>& range,
                                                 const size_t over_integrate = 0,
                                                 const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::
      make_unique<L2GlobalProlongationLocalizableOperator<typename RangeSpaceType::GridLayerType,
                                                          ConstDiscreteFunction<SourceSpaceType, SourceVectorType>,
                                                          DiscreteFunction<RangeSpaceType, RangeVectorType>>>(
          over_integrate, range.space().grid_layer(), source, range, param);
} // ... make_global_l2_prolongation_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2GlobalProlongationOperator
    : public OperatorInterface<internal::L2GlobalProlongationOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2GlobalProlongationOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2GlobalProlongationOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2GlobalProlongationOperator(const size_t over_integrate, GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  L2GlobalProlongationOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
  {
  }

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source,
             DiscreteFunction<RS, RV>& range,
             const XT::Common::Parameter& param = {}) const
  {
    L2GlobalProlongationLocalizableOperator<GridLayerType, ConstDiscreteFunction<SS, SV>, DiscreteFunction<RS, RV>> op(
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
}; // class L2GlobalProlongationOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2GlobalProlongationOperator<GridLayerType>>>::type
make_global_l2_prolongation_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2GlobalProlongationOperator<GridLayerType>>(over_integrate, grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_GLOBAL_HH
