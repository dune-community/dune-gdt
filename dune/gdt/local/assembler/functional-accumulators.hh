// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ACCUMULATORS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ACCUMULATORS_HH

#include <memory>

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/functionals/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \note See also LocalElementFunctionalInterface for a description of some of the template arguments.
 *
 * \sa make_local_element_functional_accumulator
 * \sa LocalElementFunctionalInterface
 * \sa LocalizableFunctionalBase
 */
template <class GridView, size_t r = 1, size_t rC = 1, class R = double, class F = double>
class LocalElementFunctionalAccumulator
  : public XT::Grid::ElementFunctor<GridView>
  , public XT::Common::ThreadResultPropagator<LocalElementFunctionalAccumulator<GridView, r, rC, R, F>, F>
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalElementFunctionalAccumulator;
  using BaseType = XT::Grid::ElementFunctor<GridView>;
  using Propagator = XT::Common::ThreadResultPropagator<LocalElementFunctionalAccumulator<GridView, r, rC, R, F>, F>;
  friend Propagator;

public:
  using E = XT::Grid::extract_entity_t<GridView>;
  using typename BaseType::ElementType;
  using ResultType = F;

  using SourceType = XT::Functions::GridFunctionInterface<E, r, rC, R>;
  using LocalFunctionalType = LocalElementFunctionalInterface<E, r, rC, R, F>;

  LocalElementFunctionalAccumulator(const LocalFunctionalType& local_functional,
                                    const SourceType& source,
                                    const XT::Common::Parameter& param = {})
    : BaseType()
    , Propagator(this)
    , local_functional_(local_functional.copy())
    , source_(source)
    , result_(0)
    , param_(param)
    , local_source_(source_.local_function())
  {}

  LocalElementFunctionalAccumulator(const ThisType& other)
    : BaseType(other)
    , Propagator(other)
    , local_functional_(other.local_functional_->copy())
    , source_(other.source_)
    , result_(0)
    , param_(other.param_)
    , local_source_(source_.local_function())
  {}

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

  void apply_local(const ElementType& element) override final
  {
    local_source_->bind(element);
    DUNE_THROW_IF(
        local_source_->size() != 1, Exceptions::assembler_error, "local_source_->size() = " << local_source_->size());
    local_functional_->apply(*local_source_, functional_value_, param_);
    result_ += functional_value_[0];
  }

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
  const std::unique_ptr<LocalFunctionalType> local_functional_;
  const SourceType& source_;
  ResultType result_;
  const XT::Common::Parameter param_;
  std::unique_ptr<typename SourceType::LocalFunctionType> local_source_;
  DynamicVector<F> functional_value_;
}; // class LocalElementFunctionalAccumulator


/**
 * \sa LocalElementFunctionalAccumulator
 */
template <class GridView, class E, size_t r, size_t rC, class R, class F>
std::enable_if_t<XT::Grid::is_view<GridView>::value && std::is_same<E, XT::Grid::extract_entity_t<GridView>>::value,
                 std::unique_ptr<LocalElementFunctionalAccumulator<GridView, r, rC, R, F>>>
make_local_element_functional_accumulator(const LocalElementFunctionalInterface<E, r, rC, R, F>& local_functional,
                                          const XT::Functions::GridFunctionInterface<E, r, rC, R>& source,
                                          const XT::Common::Parameter& param = {})
{
  return std::make_unique<LocalElementFunctionalAccumulator<GridView, r, rC, R, F>>(local_functional, source, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ACCUMULATORS_HH
