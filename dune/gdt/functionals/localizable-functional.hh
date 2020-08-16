// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_FUNCTIONALS_LOCALIZABLE_FUNCTIONAL_HH
#define DUNE_GDT_FUNCTIONALS_LOCALIZABLE_FUNCTIONAL_HH

#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/local/assembler/functional-accumulators.hh>
#include <dune/gdt/local/functionals/interfaces.hh>

namespace Dune {
namespace GDT {


template <class GridView,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double,
          class Field = double>
class LocalizableFunctionalBase : public XT::Grid::Walker<GridView>
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalizableFunctionalBase;
  using BaseType = XT::Grid::Walker<GridView>;

public:
  using GV = GridView;
  using GridViewType = GridView;
  using E = XT::Grid::extract_entity_t<GV>;

  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;
  using R = RangeField;
  using F = Field;
  using SourceType = XT::Functions::GridFunctionInterface<E, r, rC, R>;

  using typename BaseType::I;
  using ElementFilterType = XT::Grid::ElementFilter<GridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<GridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<GridViewType>;
  using ApplyOnAllIntersection = XT::Grid::ApplyOn::AllIntersections<GridViewType>;

  LocalizableFunctionalBase(GridViewType assembly_grid_view, const SourceType& src)
    : BaseType(assembly_grid_view)
    , source_(src.copy_as_grid_function())
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  const SourceType& source() const
  {
    return *source_;
  }

  using BaseType::append;

  ThisType& append(const LocalElementFunctionalInterface<E, r, rC, R, F>& local_functional,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    local_accumulators_.emplace_back(make_local_element_functional_accumulator<GV>(local_functional, &source_, param));
    this->append(*local_accumulators_.back(), filter);
    return *this;
  }

  void assemble(const bool use_tbb = false)
  {
    if (assembled_)
      return;
    // This clears all appended functionals, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

  Field result() const
  {
    Field ret = 0;
    for (const auto& local_accumulators : local_accumulators_)
      ret += local_accumulators->result();
    return ret;
  }

protected:
  const std::unique_ptr<SourceType> source_;
  bool assembled_;
  std::list<std::unique_ptr<LocalElementFunctionalAccumulator<GridView, r, rC, R, F>>> local_accumulators_;
}; // class LocalizableFunctionalBase


template <class GV, size_t r, size_t rC, class R>
std::enable_if_t<XT::Grid::is_view<GV>::value, LocalizableFunctionalBase<GV, r, rC, R>> make_localizable_functional(
    GV assembly_grid_view, const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, R>& source)
{
  return LocalizableFunctionalBase<GV, r, rC, R>(assembly_grid_view, source);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_LOCALIZABLE_FUNCTIONAL_HH
