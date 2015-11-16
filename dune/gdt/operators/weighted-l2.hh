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


// //////////////////////// //
// WeightedL2MatrixOperator //
// //////////////////////// //

template <class WeightFunctionType, class RangeSpace,
          class Matrix   = typename Stuff::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridView = typename RangeSpace::GridViewType, class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class WeightedL2MatrixOperator
    : public MatrixOperatorDefault<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume>
{
  typedef MatrixOperatorDefault<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume> BaseType;
  typedef LocalVolumeIntegralOperator<LocalEvaluation::Product<WeightFunctionType>> LocalWeightedL2OperatorType;

public:
  template <class... Args>
  explicit WeightedL2MatrixOperator(const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(weight)
  {
    this->add(local_weighted_l2_operator_);
  }

  template <class... Args>
  explicit WeightedL2MatrixOperator(const size_t over_integrate, const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(over_integrate, weight)
  {
    this->add(local_weighted_l2_operator_);
  }

private:
  const LocalWeightedL2OperatorType local_weighted_l2_operator_;
}; // class WeightedL2MatrixOperator


// //////////////////////////////// //
// make_weighted_l2_matrix_operator //
// //////////////////////////////// //

// without matrix

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space, grid_view);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType, class GridViewType>
typename std::
    enable_if<Stuff::LA::is_matrix<MatrixType>::value && Stuff::is_localizable_function<WeightFunctionType>::value
                  && is_space<SpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
      over_integrate, weight, space, grid_view);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, range_space, source_space, grid_view);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value
                            && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space, const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
                                                   RangeSpaceType,
                                                   MatrixType,
                                                   GridViewType,
                                                   SourceSpaceType>>(
      over_integrate, weight, range_space, source_space, grid_view);
}

// with matrix

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space, grid_view of the space is
 *        used).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, matrix, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType, class GridViewType>
typename std::
    enable_if<Stuff::is_localizable_function<WeightFunctionType>::value && Stuff::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
      over_integrate, weight, matrix, space, grid_view);
}

/**
 * \brief Creates a weighted L2 matrix operator.
 */
template <class WeightFunctionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix,
                                 const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                 const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
                                                   RangeSpaceType,
                                                   MatrixType,
                                                   GridViewType,
                                                   SourceSpaceType>>(
      over_integrate, weight, matrix, range_space, source_space, grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
