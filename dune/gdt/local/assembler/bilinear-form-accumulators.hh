// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_BILINEAR_FORM_ACCUMULATORS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_BILINEAR_FORM_ACCUMULATORS_HH

#include <memory>

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \note See also LocalElementBilinearFormInterface for a description of some of the template arguments.
 *
 * \sa make_local_element_bilinear_form_accumulator
 * \sa LocalElementBilinearFormInterface
 * \sa LocalizableBilinearFormBase
 */
template <class GV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class R = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = SF>
class LocalElementBilinearFormAccumulator
  : public XT::Grid::ElementFunctor<GV>
  , public XT::Common::ThreadResultPropagator<LocalElementBilinearFormAccumulator<GV, s_r, s_rC, SF, R, r_r, r_rC, RF>,
                                              R>
{
  static_assert(XT::Grid::is_view<GV>::value, "");

  using ThisType = LocalElementBilinearFormAccumulator<GV, s_r, s_rC, SF, R, r_r, r_rC, RF>;
  using BaseType = XT::Grid::ElementFunctor<GV>;
  using Propagator =
      XT::Common::ThreadResultPropagator<LocalElementBilinearFormAccumulator<GV, s_r, s_rC, SF, R, r_r, r_rC, RF>, R>;
  friend Propagator;

public:
  using E = XT::Grid::extract_entity_t<GV>;
  using typename BaseType::ElementType;
  using ResultType = R;

  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SF>;
  using RangeType = XT::Functions::GridFunctionInterface<E, r_r, r_rC, RF>;
  using LocalBilinearFormType = LocalElementBilinearFormInterface<E, s_r, s_rC, SF, R, r_r, r_rC, RF>;

  LocalElementBilinearFormAccumulator(const LocalBilinearFormType& local_bilinear_form,
                                      const SourceType& source,
                                      const RangeType& range,
                                      const XT::Common::Parameter& param = {})
    : BaseType()
    , Propagator(this)
    , local_bilinear_form_(local_bilinear_form.copy())
    , source_(source)
    , range_(range)
    , result_(0)
    , param_(param)
    , local_source_(source_.local_function())
    , local_range_(range_.local_function())
  {}

  LocalElementBilinearFormAccumulator(const ThisType& other)
    : BaseType(other)
    , Propagator(other)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
    , source_(other.source_)
    , range_(other.range_)
    , result_(0)
    , param_(other.param_)
    , local_source_(source_.local_function())
    , local_range_(range_.local_function())
  {}

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

  void apply_local(const ElementType& element) override final
  {
    local_source_->bind(element);
    local_range_->bind(element);
    DUNE_THROW_IF(
        local_source_->size() != 1, Exceptions::assembler_error, "local_source_->size() = " << local_source_->size());
    DUNE_THROW_IF(
        local_range_->size() != 1, Exceptions::assembler_error, "local_range_->size() = " << local_range_->size());
    local_bilinear_form_->apply2(*local_source_, *local_range_, bilinear_form_value_, param_);
    result_ += bilinear_form_value_[0][0];
  } // ... apply_local(...)

  void finalize() override final
  {
    Propagator::finalize_imp();
  }

  const ResultType& result() const
  {
    return result_;
  }

protected:
  void set_result(const ResultType& res)
  {
    result_ = res;
  }

private:
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  const SourceType& source_;
  const RangeType& range_;
  ResultType result_;
  const XT::Common::Parameter param_;
  std::unique_ptr<typename SourceType::LocalFunctionType> local_source_;
  std::unique_ptr<typename RangeType::LocalFunctionType> local_range_;
  DynamicMatrix<R> bilinear_form_value_;
}; // class LocalElementBilinearFormAccumulator


/**
 * \sa LocalElementBilinearFormAccumulator
 */
template <class GridView, class E, size_t s_r, size_t s_rC, class SF, size_t r_r, size_t r_rC, class RF, class R>
std::enable_if_t<XT::Grid::is_view<GridView>::value && std::is_same<E, XT::Grid::extract_entity_t<GridView>>::value,
                 std::unique_ptr<LocalElementBilinearFormAccumulator<GridView, s_r, s_rC, SF, R, r_r, r_rC, RF>>>
make_local_element_bilinear_form_accumulator(
    const LocalElementBilinearFormInterface<E, s_r, s_rC, SF, R, r_r, r_rC, RF>& local_bilinear_form,
    const XT::Functions::GridFunctionInterface<E, s_r, s_rC, SF>& source,
    const XT::Functions::GridFunctionInterface<E, r_r, r_rC, RF>& range,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<LocalElementBilinearFormAccumulator<GridView, s_r, s_rC, SF, R, r_r, r_rC, RF>>(
      local_bilinear_form, source, range, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_BILINEAR_FORM_ACCUMULATORS_HH
