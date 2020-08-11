// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_LOCAL_OPERATORS_INDICATOR_HH
#define DUNE_GDT_LOCAL_OPERATORS_INDICATOR_HH

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class SV,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class RF = SF,
          class RGV = SGV,
          class RV = SV>
class LocalElementBilinearFormIndicatorOperator
  : public LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, 1, 1, RF, RGV, RV>
{
public:
  using ThisType = LocalElementBilinearFormIndicatorOperator;
  using BaseType = LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, 1, 1, RF, RGV, RV>;

  using typename BaseType::E;
  using typename BaseType::LocalRangeType;
  using typename BaseType::SourceType;

  using LocalBilinearFormType = LocalElementBilinearFormInterface<E, s_r, s_rC, SF, RF, s_r, s_rC, SF>;

  // When using this constructor, source has to be set by a call to with_source before calling apply
  LocalElementBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form)
    : BaseType(1)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalElementBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form, const SourceType& src)
    : BaseType(src)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalElementBilinearFormIndicatorOperator(const ThisType& other)
    : BaseType(other)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const override final
  {
    const auto& local_source = *(this->local_sources_[0]);
    DUNE_THROW_IF(local_source.size(param) != 1,
                  Exceptions::operator_error,
                  "u->size(" << param << ") = " << local_source.size(param));
    this->local_bilinear_form_->apply2(local_source, local_source, local_bilinearform_result_, param);
    if (local_range.space().type() == SpaceType::finite_volume)
      local_range.dofs()[0] = local_bilinearform_result_[0][0];
    else {
      const auto& value = local_bilinearform_result_[0][0];
      local_range.basis().interpolate([&value](const auto&) { return value; }, 0, local_range.dofs());
    }
  } // ... apply(...)

private:
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  mutable DynamicMatrix<RF> local_bilinearform_result_;
}; // class LocalElementBilinearFormIndicatorOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INDICATOR_HH
