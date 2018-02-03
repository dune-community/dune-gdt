// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_LOCAL_ACCUMULATROS_HH
#define DUNE_GDT_ASSEMBLER_LOCAL_ACCUMULATROS_HH

#include <dune/xt/grid/walker/wrapper.hh>

#include <dune/gdt/local/operators/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \todo \attention Rename LocalVolumeTwoFormAccumulatorFunctor -> LocalVolumeTwoFormAccumulator once the latter is
 *                  removed!
 */
template <class GL,
          class TestFunctionType,
          class AnsatzFunctionType,
          class FieldType = typename TestFunctionType::RangeFieldType>
class LocalVolumeTwoFormAccumulatorFunctor : public XT::Grid::internal::Codim0ReturnObject<GL, FieldType>
{
  static_assert(XT::Functions::is_localizable_function<TestFunctionType>::value, "");
  static_assert(XT::Functions::is_localizable_function<AnsatzFunctionType>::value, "");

  typedef LocalVolumeTwoFormAccumulatorFunctor<GL, TestFunctionType, AnsatzFunctionType, FieldType> ThisType;
  typedef XT::Grid::internal::Codim0ReturnObject<GL, FieldType> BaseType;

public:
  typedef LocalVolumeTwoFormInterface<typename TestFunctionType::LocalfunctionType,
                                      typename AnsatzFunctionType::LocalfunctionType,
                                      FieldType>
      LocalVolumeTwoFormType;
  typedef typename BaseType::GridLayerType GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  LocalVolumeTwoFormAccumulatorFunctor(const GridLayerType& grid_layer,
                                       const LocalVolumeTwoFormType& local_volume_two_form,
                                       const TestFunctionType& test_function,
                                       const AnsatzFunctionType& ansatz_function,
                                       FieldType& res,
                                       const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where)
    : grid_layer_(grid_layer)
    , local_volume_two_form_(local_volume_two_form)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , final_result_(res)
    , where_(where)
    , result_per_thread_(0.)
    , finalized_(false)
  {
  }

  void prepare() override final
  {
    final_result_ = 0.;
    *result_per_thread_ = 0.;
  }

  bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return where_->apply_on(grid_layer, entity);
  }

  FieldType compute_locally(const EntityType& entity) override final
  {
    DynamicMatrix<FieldType> local_twoform_result(1, 1, 0.); // \todo: make mutable member, after SMP refactor
    this->local_volume_two_form_.apply2(
        *test_function_.local_function(entity), *ansatz_function_.local_function(entity), local_twoform_result);
    return local_twoform_result[0][0];
  } // ... compute_locally(...)

  void apply_local(const EntityType& entity) override final
  {
    *result_per_thread_ += compute_locally(entity);
  }

  void finalize() override final
  {
    if (!finalized_) {
      final_result_ = result_per_thread_.sum();
      final_result_ = grid_layer_.comm().sum(final_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  FieldType result() const override final
  {
    DUNE_THROW_IF(!finalized_, XT::Common::Exceptions::you_are_using_this_wrong, "Call finalize() first!");
    return final_result_;
  }

private:
  const GridLayerType& grid_layer_;
  const LocalVolumeTwoFormType& local_volume_two_form_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  FieldType& final_result_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  Dune::XT::Common::PerThreadValue<FieldType> result_per_thread_;
  bool finalized_;
}; // class LocalVolumeTwoFormAccumulatorFunctor


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_LOCAL_ACCUMULATROS_HH
