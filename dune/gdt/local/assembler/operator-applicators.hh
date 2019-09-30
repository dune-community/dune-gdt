// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_APPLICATORS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_APPLICATORS_HH

#include <memory>

#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/operators/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \note See also LocalElementOperatorInterface for a description of the template arguments.
 *
 * \sa make_local_element_operator_applicator
 * \sa LocalElementOperatorInterface
 * \sa LocalizableOperatorBase
 */
template <class AssemblyGridView,
          class SV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = SF,
          class RGV = SGV,
          class RV = SV>
class LocalElementOperatorApplicator : public XT::Grid::ElementFunctor<AssemblyGridView>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_vector<SV>::value, "");
  static_assert(XT::Grid::is_view<SGV>::value, "");
  static_assert(XT::Grid::is_view<RGV>::value, "");
  static_assert(XT::LA::is_vector<RV>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<AssemblyGridView>, XT::Grid::extract_entity_t<SGV>>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<AssemblyGridView>, XT::Grid::extract_entity_t<RGV>>::value, "");

  using ThisType = LocalElementOperatorApplicator;
  using BaseType = XT::Grid::ElementFunctor<AssemblyGridView>;

public:
  using typename BaseType::ElementType;

  using SourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>;
  using RangeType = DiscreteFunction<RV, RGV, r_r, r_rC, RF>;
  using LocalOperatorType = LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>;

  LocalElementOperatorApplicator(const LocalOperatorType& local_operator,
                                 RangeType& range,
                                 const XT::Common::Parameter& param = {})
    : local_operator_(local_operator.copy())
    , range_(range)
    , param_(param)
    , local_range_(range_.local_discrete_function())
  {}

  LocalElementOperatorApplicator(const ThisType& other)
    : LocalElementOperatorApplicator(*other.local_operator_, other.range_, other.param_)
  {}

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const ElementType& element) override final
  {
    local_range_->bind(element);
    local_operator_->apply(*local_range_, param_);
  }

private:
  const std::unique_ptr<LocalOperatorType> local_operator_;
  RangeType& range_;
  const XT::Common::Parameter param_;
  std::unique_ptr<typename RangeType::LocalDiscreteFunctionType> local_range_;
}; // class LocalElementOperatorApplicator


/**
 * \sa LocalElementOperatorApplicator
 */
template <class SV, class GV, size_t s_r, size_t s_rC, class SF, size_t r_r, size_t r_rC, class RF, class RV>
std::unique_ptr<LocalElementOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, RV>>
make_local_element_operator_applicator(
    const LocalElementOperatorInterface<SV, GV, s_r, s_rC, SF, r_r, r_rC, RF, GV, RV>& local_operator,
    DiscreteFunction<RV, GV, r_r, r_rC, RF>& range,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<LocalElementOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, RV>>(
      local_operator, range, param);
}

/**
 * \sa LocalElementOperatorApplicator
 */
template <class AssemblyGridView,
          class SV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV>
std::unique_ptr<LocalElementOperatorApplicator<AssemblyGridView, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>
make_local_element_operator_applicator(
    const LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>& local_operator,
    DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<
      LocalElementOperatorApplicator<AssemblyGridView, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>(
      local_operator, range, param);
}


/**
 * \note See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa make_local_intersection_operator_applicator
 * \sa LocalIntersectionOperatorInterface
 * \sa LocalizableOperatorBase
 */
template <class AssemblyGridView,
          class SV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RF = SF,
          class IRGV = SGV,
          class IRV = SV,
          class ORGV = SGV,
          class ORV = SV>
class LocalIntersectionOperatorApplicator : public XT::Grid::IntersectionFunctor<AssemblyGridView>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_vector<SV>::value, "");
  static_assert(XT::Grid::is_view<SGV>::value, "");
  static_assert(XT::Grid::is_view<IRGV>::value, "");
  static_assert(XT::LA::is_vector<IRV>::value, "");
  static_assert(XT::Grid::is_view<ORGV>::value, "");
  static_assert(XT::LA::is_vector<ORV>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<AssemblyGridView>, XT::Grid::extract_entity_t<SGV>>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<AssemblyGridView>, XT::Grid::extract_entity_t<IRGV>>::value,
                "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<AssemblyGridView>, XT::Grid::extract_entity_t<ORGV>>::value,
                "");

  using ThisType = LocalIntersectionOperatorApplicator<AssemblyGridView,
                                                       SV,
                                                       s_r,
                                                       s_rC,
                                                       SF,
                                                       SGV,
                                                       r_r,
                                                       r_rC,
                                                       RF,
                                                       IRGV,
                                                       IRV,
                                                       ORGV,
                                                       ORV>;
  using BaseType = XT::Grid::IntersectionFunctor<AssemblyGridView>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;

  using SourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>;
  using InsideRangeType = DiscreteFunction<IRV, IRGV, r_r, r_rC, RF>;
  using OutsideRangeType = DiscreteFunction<ORV, ORGV, r_r, r_rC, RF>;
  using LocalOperatorType =
      LocalIntersectionOperatorInterface<IntersectionType, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, IRGV, IRV, ORGV, ORV>;

  LocalIntersectionOperatorApplicator(const LocalOperatorType& local_operator,
                                      InsideRangeType& range_inside,
                                      OutsideRangeType& range_outside,
                                      const XT::Common::Parameter& param = {})
    : local_operator_(local_operator.copy())
    , range_inside_(range_inside)
    , range_outside_(range_outside)
    , param_(param)
    , local_range_inside_(range_inside_.local_discrete_function())
    , local_range_outside_(range_outside_.local_discrete_function())
  {}

  LocalIntersectionOperatorApplicator(const ThisType& other)
    : LocalIntersectionOperatorApplicator(
          *other.local_operator_, other.range_inside_, other.range_outside_, other.param_)
  {}

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    local_range_inside_->bind(inside_element);
    local_range_outside_->bind(outside_element);
    local_operator_->bind(intersection);
    local_operator_->apply(intersection, *local_range_inside_, *local_range_outside_, param_);
  }

private:
  const std::unique_ptr<LocalOperatorType> local_operator_;
  InsideRangeType& range_inside_;
  OutsideRangeType& range_outside_;
  const XT::Common::Parameter param_;
  std::unique_ptr<typename InsideRangeType::LocalDiscreteFunctionType> local_range_inside_;
  std::unique_ptr<typename OutsideRangeType::LocalDiscreteFunctionType> local_range_outside_;
}; // class LocalIntersectionOperatorApplicator


/**
 * \{
 * \name Variants of make_local_intersection_operator_applicator without a manually specified AssemblyGridView.
 */

template <class I,
          class SV,
          class GV,
          size_t s_r,
          size_t s_rC,
          class SF,
          size_t r_r,
          size_t r_rC,
          class RF,
          class IRV,
          class ORV>
std::enable_if_t<
    std::is_same<XT::Grid::extract_intersection_t<GV>, I>::value,
    std::unique_ptr<LocalIntersectionOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, IRV, GV, ORV>>>
make_local_intersection_operator_applicator(
    const LocalIntersectionOperatorInterface<I, SV, GV, s_r, s_rC, SF, r_r, r_rC, RF, GV, IRV, GV, ORV>& local_operator,
    DiscreteFunction<IRV, GV, r_r, r_rC, RF>& range_inside,
    DiscreteFunction<ORV, GV, r_r, r_rC, RF>& range_outside,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<
      LocalIntersectionOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, IRV, GV, ORV>>(
      local_operator, range_inside, range_outside, param);
} // ... make_local_intersection_operator_applicator(...)

template <class I, class SV, class GV, size_t s_r, size_t s_rC, class SF, size_t r_r, size_t r_rC, class RF, class RV>
std::enable_if_t<std::is_same<XT::Grid::extract_intersection_t<GV>, I>::value,
                 std::unique_ptr<LocalIntersectionOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, RV>>>
make_local_intersection_operator_applicator(
    const LocalIntersectionOperatorInterface<I, SV, GV, s_r, s_rC, SF, r_r, r_rC, RF, GV, RV>& local_operator,
    DiscreteFunction<RV, GV, r_r, r_rC, RF>& range,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<LocalIntersectionOperatorApplicator<GV, SV, s_r, s_rC, SF, GV, r_r, r_rC, RF, GV, RV>>(
      local_operator, range, range, param);
}

/**
 * \{
 * \name Variants of make_local_intersection_operator_applicator with a manually specified AssemblyGridView.
 */

template <class AssemblyGridView,
          class I,
          class SV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          size_t r_r,
          size_t r_rC,
          class RF,
          class IRGV,
          class IRV,
          class ORGV,
          class ORV>
std::enable_if_t<std::is_same<XT::Grid::extract_intersection_t<AssemblyGridView>, I>::value,
                 std::unique_ptr<LocalIntersectionOperatorApplicator<AssemblyGridView,
                                                                     SV,
                                                                     s_r,
                                                                     s_rC,
                                                                     SF,
                                                                     SGV,
                                                                     r_r,
                                                                     r_rC,
                                                                     RF,
                                                                     IRGV,
                                                                     IRV,
                                                                     ORGV,
                                                                     ORV>>>
make_local_intersection_operator_applicator(
    const LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, IRGV, IRV, ORGV, ORV>&
        local_operator,
    DiscreteFunction<IRV, IRGV, r_r, r_rC, RF>& range_inside,
    DiscreteFunction<ORV, ORGV, r_r, r_rC, RF>& range_outside,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<LocalIntersectionOperatorApplicator<AssemblyGridView,
                                                              SV,
                                                              s_r,
                                                              s_rC,
                                                              SF,
                                                              SGV,
                                                              r_r,
                                                              r_rC,
                                                              RF,
                                                              IRGV,
                                                              IRV,
                                                              ORGV,
                                                              ORV>>(local_operator, range_inside, range_outside, param);
} // ... make_local_intersection_operator_applicator(...)

template <class AssemblyGridView,
          class I,
          class SV,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV>
std::enable_if_t<
    std::is_same<XT::Grid::extract_intersection_t<AssemblyGridView>, I>::value,
    std::unique_ptr<
        LocalIntersectionOperatorApplicator<AssemblyGridView, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>>
make_local_intersection_operator_applicator(
    const LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>& local_operator,
    DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
    const XT::Common::Parameter& param = {})
{
  return std::make_unique<
      LocalIntersectionOperatorApplicator<AssemblyGridView, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>(
      local_operator, range, range, param);
}

/// \}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_APPLICATORS_HH
