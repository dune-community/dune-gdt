// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_OPERATORS_WEIGHTED_L2_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localoperator/integrals.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


// //////////////////////////// //
// WeightedL2LocalizableProduct //
// //////////////////////////// //

template <class WeightFunctionType, class GridView, class Range, class Source = Range,
          class Field                                                         = typename Range::RangeFieldType>
class WeightedL2LocalizableProduct : public LocalizableProductDefault<GridView, Range, Source, Field>
{
  typedef LocalizableProductDefault<GridView, Range, Source, Field> BaseType;
  typedef LocalVolumeIntegralOperator<LocalEvaluation::Product<WeightFunctionType>> LocalWeightedL2OperatorType;

public:
  template <class... Args>
  WeightedL2LocalizableProduct(const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(weight)
  {
    this->add(local_weighted_l2_operator_);
  }

  template <class... Args>
  WeightedL2LocalizableProduct(const size_t over_integrate, const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(over_integrate, weight)
  {
    this->add(local_weighted_l2_operator_);
  }

private:
  const LocalWeightedL2OperatorType local_weighted_l2_operator_;
}; // class WeightedL2LocalizableProduct


// //////////////////////////////////// //
// make_weighted_l2_localizable_product //
// //////////////////////////////////// //

/**
 * \sa WeightedL2LocalizableProduct
 */
template <class WeightFunctionType, class GridViewType, class RangeType, class SourceType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value
                            && Stuff::is_localizable_function<RangeType>::value
                            && Stuff::is_localizable_function<SourceType>::value,
                        std::unique_ptr<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType,
                                                                     SourceType>>>::type
make_weighted_l2_localizable_product(const WeightFunctionType& weight, const GridViewType& grid_view,
                                     const RangeType& range, const SourceType& source, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType, SourceType>>(
      over_integrate, weight, grid_view, range, source);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
