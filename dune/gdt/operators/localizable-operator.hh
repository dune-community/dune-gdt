// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
#define DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH

#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/assembler/operator-applicators.hh>
#include <dune/gdt/local/operators/interfaces.hh>

namespace Dune {
namespace GDT {


template <class AssemblyGridView,
          class SourceVector,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceRangeField = double,
          class SourceGridView = AssemblyGridView,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeRangeField = SourceRangeField,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
class LocalizableOperatorBase : public XT::Grid::Walker<AssemblyGridView>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_vector<SourceVector>::value, "");
  static_assert(XT::Grid::is_view<SourceGridView>::value, "");
  static_assert(XT::Grid::is_view<RangeGridView>::value, "");
  static_assert(XT::LA::is_vector<RangeVector>::value, "");

  using ThisType = LocalizableOperatorBase<AssemblyGridView,
                                           SourceVector,
                                           source_range_dim,
                                           source_range_dim_cols,
                                           SourceRangeField,
                                           SourceGridView,
                                           range_range_dim,
                                           range_range_dim_cols,
                                           RangeRangeField,
                                           RangeGridView,
                                           RangeVector>;
  using BaseType = XT::Grid::Walker<AssemblyGridView>;

public:
  using AGV = AssemblyGridView;
  using AssemblyGridViewType = AssemblyGridView;

  using SV = SourceVector;
  using SGV = SourceGridView;
  static const constexpr size_t s_r = source_range_dim;
  static const constexpr size_t s_rC = source_range_dim_cols;
  using SR = SourceRangeField;
  using SourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SR>;

  using RV = RangeVector;
  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_range_dim;
  static const constexpr size_t r_rC = range_range_dim_cols;
  using RR = RangeRangeField;
  using RangeType = DiscreteFunction<RV, RGV, r_r, r_rC, RR>;

  using typename BaseType::I;
  using ElementFilterType = XT::Grid::ElementFilter<AssemblyGridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AssemblyGridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AssemblyGridViewType>;
  using ApplyOnAllIntersection = XT::Grid::ApplyOn::AllIntersections<AssemblyGridViewType>;

  LocalizableOperatorBase(AssemblyGridViewType assembly_grid_view, const SourceType& src, RangeType& rng)
    : BaseType(assembly_grid_view)
    , source_(src)
    , range_(rng)
    , assembled_(false)
  {
    // to detect assembly
    this->append([&](const auto&) { assembled_ = true; });
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

  using BaseType::append;

  ThisType& append(const LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SR, r_r, r_rC, RR, RGV, RV>& local_operator,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    this->append(make_local_element_operator_applicator(local_operator, source_, range_, param).release(), filter);
    return *this;
  }

  ThisType&
  append(const LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SR, r_r, r_rC, RR, RGV, RV>& local_operator,
         const XT::Common::Parameter& param = {},
         const IntersectionFilterType& filter = ApplyOnAllIntersection())
  {
    this->append(make_local_intersection_operator_applicator(local_operator, source_, range_, param).release(), filter);
    return *this;
  }

  void assemble(const bool use_tbb = false)
  {
    if (assembled_)
      return;
    // This clears all appended operators, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

protected:
  const SourceType& source_;
  RangeType& range_;
  bool assembled_;
}; // class LocalizableOperatorBase


template <class AGV,
          class SV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class SGV,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV>
std::enable_if_t<XT::Grid::is_layer<AGV>::value,
                 LocalizableOperatorBase<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>
make_localizable_operator(AGV assembly_grid_view,
                          const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                          DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range)
{
  return LocalizableOperatorBase<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>(
      assembly_grid_view, source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
