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


template <class I,
          class SV,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class RF = SF,
          class IRGV = SGV,
          class IRV = SV,
          class ORGV = SGV,
          class ORV = SV>
class LocalCouplingIntersectionBilinearFormIndicatorOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, 1, 1, RF, IRGV, IRV, ORGV, ORV>
{
public:
  using ThisType = LocalCouplingIntersectionBilinearFormIndicatorOperator;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, 1, 1, RF, IRGV, IRV, ORGV, ORV>;

  using typename BaseType::E;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceType;

  using LocalBilinearFormType = LocalCouplingIntersectionBilinearFormInterface<I, s_r, s_rC, SF, RF, s_r, s_rC, SF>;

  // When using this constructor, source has to be set by a call to with_source before calling apply
  LocalCouplingIntersectionBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form)
    : BaseType(2)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalCouplingIntersectionBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form,
                                                         const SourceType& src)
    : BaseType(src, 2)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalCouplingIntersectionBilinearFormIndicatorOperator(const ThisType& other)
    : BaseType(other)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& local_source_in = *(this->local_sources_[0]);
    const auto& local_source_out = *(this->local_sources_[1]);
    DUNE_THROW_IF(local_source_in.size(param) != 1,
                  Exceptions::operator_error,
                  "u_in.size(" << param << ") = " << local_source_in.size(param));
    DUNE_THROW_IF(local_source_out.size(param) != 1,
                  Exceptions::operator_error,
                  "u_out.size(" << param << ") = " << local_source_out.size(param));
    this->local_bilinear_form_->apply2(this->intersection(),
                                       local_source_in,
                                       local_source_in,
                                       local_source_out,
                                       local_source_out,
                                       local_bilinearform_result_in_in_,
                                       local_bilinearform_result_in_out_,
                                       local_bilinearform_result_out_in_,
                                       local_bilinearform_result_out_out_,
                                       param);
    DUNE_THROW_IF(local_range_inside.space().type() != SpaceType::finite_volume_skeleton,
                  Exceptions::operator_error,
                  "Not implemented yet for non-skeleton spaces!");
    // For the skeleton FV space, we know how to get the global index of this intersection:
    const auto& element = local_source_in.element();
    const auto& intersection_mapper = local_range_inside.space().mapper();
    const auto intersection_index = intersection_mapper.global_index(element, this->intersection().indexInInside());
    local_range_inside.dofs().global()[intersection_index] =
        local_bilinearform_result_in_in_[0][0] + local_bilinearform_result_in_out_[0][0]
        + local_bilinearform_result_out_in_[0][0] + local_bilinearform_result_out_out_[0][0];
  } // ... apply(...)

private:
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  mutable DynamicMatrix<RF> local_bilinearform_result_in_in_;
  mutable DynamicMatrix<RF> local_bilinearform_result_in_out_;
  mutable DynamicMatrix<RF> local_bilinearform_result_out_in_;
  mutable DynamicMatrix<RF> local_bilinearform_result_out_out_;
}; // class LocalCouplingIntersectionBilinearFormIndicatorOperator


template <class I,
          class SV,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class RF = SF,
          class IRGV = SGV,
          class IRV = SV,
          class ORGV = SGV,
          class ORV = SV>
class LocalIntersectionBilinearFormIndicatorOperator
  : public LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, 1, 1, RF, IRGV, IRV, ORGV, ORV>
{
public:
  using ThisType = LocalIntersectionBilinearFormIndicatorOperator;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, 1, 1, RF, IRGV, IRV, ORGV, ORV>;

  using typename BaseType::E;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;
  using typename BaseType::SourceType;

  using LocalBilinearFormType = LocalIntersectionBilinearFormInterface<I, s_r, s_rC, SF, RF, s_r, s_rC, SF>;

  // When using this constructor, source has to be set by a call to with_source before calling apply
  LocalIntersectionBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form)
    : BaseType(1)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalIntersectionBilinearFormIndicatorOperator(const LocalBilinearFormType& bilinear_form, const SourceType& src)
    : BaseType(src, 1)
    , local_bilinear_form_(bilinear_form.copy())
  {}

  LocalIntersectionBilinearFormIndicatorOperator(const ThisType& other)
    : BaseType(other)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& /*local_range_outside*/,
             const XT::Common::Parameter& param = {}) const override final
  {
    const auto& local_source = *(this->local_sources_[0]);
    DUNE_THROW_IF(local_source.size(param) != 1,
                  Exceptions::operator_error,
                  "u.size(" << param << ") = " << local_source.size(param));
    this->local_bilinear_form_->apply2(
        this->intersection(), local_source, local_source, local_bilinearform_result_, param);
    DUNE_THROW_IF(local_range_inside.space().type() != SpaceType::finite_volume_skeleton,
                  Exceptions::operator_error,
                  "Not implemented yet for non-skeleton spaces!");
    // For the skeleton FV space, we know how to get the global index of this intersection:
    const auto& element = local_source.element();
    const auto& intersection_mapper = local_range_inside.space().mapper();
    const auto intersection_index = intersection_mapper.global_index(
        element,
        local_bilinear_form_->inside() ? this->intersection().indexInInside() : this->intersection().indexInOutside());
    local_range_inside.dofs().global()[intersection_index] = local_bilinearform_result_[0][0];
  } // ... apply(...)

protected:
  void post_bind(const I& intrsctn) override final
  {
    if (local_bilinear_form_->inside())
      this->local_sources_[0]->bind(intrsctn.inside());
    else {
      DUNE_THROW_IF(!intrsctn.neighbor(),
                    Exceptions::operator_error,
                    "Intersection has no outside element, but bilinear form request bases to be bound to the outside!"
                        << "\nDid you forget to use the correct IntersectionFilter?");
      this->local_sources_[0]->bind(intrsctn.outside());
    }
  } // ... post_bind(...)

private:
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  mutable DynamicMatrix<RF> local_bilinearform_result_;
}; // class LocalIntersectionBilinearFormIndicatorOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INDICATOR_HH
