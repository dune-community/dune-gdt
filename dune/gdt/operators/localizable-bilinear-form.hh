// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_LOCALIZABLE_FUNCTIONAL_HH
#define DUNE_GDT_OPERATORS_LOCALIZABLE_FUNCTIONAL_HH

#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/local/assembler/bilinear-form-accumulators.hh>
#include <dune/gdt/local/bilinear-forms/interfaces.hh>

namespace Dune {
namespace GDT {


template <class GridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceRangeField = double,
          class Result = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeRangeField = SourceRangeField>
class LocalizableBilinearFormBase : public XT::Grid::Walker<GridView>
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalizableBilinearFormBase<GridView,
                                               source_range_dim,
                                               source_range_dim_cols,
                                               SourceRangeField,
                                               Result,
                                               range_range_dim,
                                               range_range_dim_cols,
                                               RangeRangeField>;
  using BaseType = XT::Grid::Walker<GridView>;

public:
  using GV = GridView;
  using GridViewType = GridView;
  using E = XT::Grid::extract_entity_t<GV>;
  using ResultType = Result;

  static constexpr size_t s_r = source_range_dim;
  static constexpr size_t s_rC = source_range_dim_cols;
  using SR = SourceRangeField;
  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SR>;

  static constexpr size_t r_r = range_range_dim;
  static constexpr size_t r_rC = range_range_dim_cols;
  using RR = RangeRangeField;
  using RangeType = XT::Functions::GridFunctionInterface<E, r_r, r_rC, RR>;

  using typename BaseType::I;
  using ElementFilterType = XT::Grid::ElementFilter<GridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<GridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<GridViewType>;
  using ApplyOnAllIntersection = XT::Grid::ApplyOn::AllIntersections<GridViewType>;

  LocalizableBilinearFormBase(GridViewType grid_view, const SourceType& src, const RangeType& rng)
    : BaseType(grid_view)
    , source_(src)
    , range_(rng)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  using BaseType::append;

  ThisType&
  append(const LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>& local_bilinear_form,
         const XT::Common::Parameter& param = {},
         const ElementFilterType& filter = ApplyOnAllElements())
  {
    local_accumulators_.emplace_back(
        make_local_element_bilinear_form_accumulator<GV>(local_bilinear_form, source_, range_, param));
    this->append(*local_accumulators_.back(), filter);
    return *this;
  }

  ThisType& append(
      std::unique_ptr<LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>>&& local_bilinear_form,
      const XT::Common::Parameter& param = {},
      const ElementFilterType& filter = ApplyOnAllElements())
  {
    local_forms_.emplace_back(std::move(local_bilinear_form));
    local_accumulators_.emplace_back(
        make_local_element_bilinear_form_accumulator<GV>(*local_forms_.back(), source_, range_, param));
    this->append(*local_accumulators_.back(), filter);
    return *this;
  }

  void assemble(const bool use_tbb = false)
  {
    if (assembled_)
      return;
    // This clears all appended bilinear forms, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

  ResultType result() const
  {
    ResultType ret = 0;
    for (const auto& local_accumulators : local_accumulators_)
      ret += local_accumulators->result();
    return ret;
  }

protected:
  const SourceType& source_;
  const RangeType& range_;
  bool assembled_;
  std::list<std::unique_ptr<LocalElementBilinearFormInterface<E, s_r, s_rC, SR, Result, r_r, r_rC, RR>>> local_forms_;
  std::list<std::unique_ptr<LocalElementBilinearFormAccumulator<GridView, s_r, s_rC, SR, ResultType, r_r, r_rC, RR>>>
      local_accumulators_;
}; // class LocalizableBilinearFormBase


template <class GV, size_t s_r, size_t s_rC, class SR, size_t r_r, size_t r_rC, class RR>
std::enable_if_t<XT::Grid::is_view<GV>::value, LocalizableBilinearFormBase<GV, s_r, s_rC, SR, double, r_r, r_rC, RR>>
make_localizable_bilinear_form(
    GV grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, s_r, s_rC, SR>& source,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r_r, r_rC, RR>& range)
{
  return LocalizableBilinearFormBase<GV, s_r, s_rC, SR, double, r_r, r_rC, RR>(grid_view, source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LOCALIZABLE_FUNCTIONAL_HH
