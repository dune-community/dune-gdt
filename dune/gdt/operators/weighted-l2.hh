// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_OPERATORS_WEIGHTED_L2_HH

#include <type_traits>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


// //////////////////////////// //
// WeightedL2LocalizableProduct //
// //////////////////////////// //

template <class WeightFunctionType,
          class GridLayer,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType>
class WeightedL2LocalizableProduct : public LocalizableProductBase<GridLayer, Range, Source, Field>
{
  typedef LocalizableProductBase<GridLayer, Range, Source, Field> BaseType;

public:
  WeightedL2LocalizableProduct(const WeightFunctionType& weight,
                               GridLayer grd_layr,
                               const Range& rng,
                               const Source& src,
                               const XT::Common::Parameter& param = {})
    : BaseType(grd_layr, rng, src)
    , local_weighted_l2_operator_(weight, param)
  {
    this->append(local_weighted_l2_operator_);
  }

  WeightedL2LocalizableProduct(const size_t over_integrate,
                               const WeightFunctionType& weight,
                               GridLayer grd_layr,
                               const Range& rng,
                               const Source& src,
                               const XT::Common::Parameter& param = {})
    : BaseType(grd_layr, rng, src)
    , local_weighted_l2_operator_(over_integrate, weight, param)
  {
    this->append(local_weighted_l2_operator_);
  }

private:
  const LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>,
                                    typename Range::LocalfunctionType,
                                    typename Source::LocalfunctionType,
                                    Field>
      local_weighted_l2_operator_;
}; // class WeightedL2LocalizableProduct


// //////////////////////////////////// //
// make_weighted_l2_localizable_product //
// //////////////////////////////////// //

/**
 * \sa WeightedL2LocalizableProduct
 */
template <class WeightFunctionType, class GridLayerType, class RangeType, class SourceType>
typename std::enable_if<
    XT::Functions::is_localizable_function<WeightFunctionType>::value && XT::Grid::is_layer<GridLayerType>::value
        && XT::Functions::is_localizable_function<RangeType>::value
        && XT::Functions::is_localizable_function<SourceType>::value,
    std::unique_ptr<WeightedL2LocalizableProduct<WeightFunctionType, GridLayerType, RangeType, SourceType>>>::type
make_weighted_l2_localizable_product(const WeightFunctionType& weight,
                                     const GridLayerType& grid_layer,
                                     const RangeType& range,
                                     const SourceType& source,
                                     const size_t over_integrate = 0,
                                     const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<
      WeightedL2LocalizableProduct<WeightFunctionType, GridLayerType, RangeType, SourceType>>(
      over_integrate, weight, grid_layer, range, source, param);
}


// //////////////////////// //
// WeightedL2MatrixOperator //
// //////////////////////// //

template <class WeightFunctionType,
          class RangeSpace,
          class Matrix = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridLayer = typename RangeSpace::GridLayerType,
          class SourceSpace = RangeSpace,
          class Field = typename Matrix::ScalarType>
class WeightedL2MatrixOperator
  : public MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::volume>
{
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::volume> BaseType;

public:
  template <class... Args>
  explicit WeightedL2MatrixOperator(const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(weight)
  {
    this->append(local_weighted_l2_operator_);
  }

  template <class... Args>
  explicit WeightedL2MatrixOperator(const size_t over_integrate, const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(over_integrate, weight)
  {
    this->append(local_weighted_l2_operator_);
  }

private:
  const LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>,
                                    typename RangeSpace::BaseFunctionSetType,
                                    typename SourceSpace::BaseFunctionSetType,
                                    Field>
      local_weighted_l2_operator_;
}; // class WeightedL2MatrixOperator


// //////////////////////////////// //
// make_weighted_l2_matrix_operator //
// //////////////////////////////// //

// without matrix

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<WeightFunctionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space, grid_layer);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType, class GridLayerType>
typename std::enable_if<
    XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<WeightFunctionType>::value
        && is_space<SpaceType>::value && XT::Grid::is_layer<GridLayerType>::value,
    std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridLayerType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 const SpaceType& space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridLayerType>>(
      over_integrate, weight, space, grid_layer);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, range_space, source_space, grid_layer);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<
    XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<WeightFunctionType>::value
        && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
        && XT::Grid::is_layer<GridLayerType>::value,
    std::unique_ptr<
        WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>(
      over_integrate, weight, range_space, source_space, grid_layer);
}

// with matrix

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space, grid_layer of the space is
 *        used).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<WeightFunctionType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 MatrixType& matrix,
                                 const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, matrix, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<
    XT::Functions::is_localizable_function<WeightFunctionType>::value && XT::LA::is_matrix<MatrixType>::value
        && is_space<SpaceType>::value && XT::Grid::is_layer<GridLayerType>::value,
    std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridLayerType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 MatrixType& matrix,
                                 const SpaceType& space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridLayerType>>(
      over_integrate, weight, matrix, space, grid_layer);
}

/**
 * \brief Creates a weighted L2 matrix operator.
 */
template <class WeightFunctionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<
    XT::Functions::is_localizable_function<WeightFunctionType>::value && XT::LA::is_matrix<MatrixType>::value
        && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
        && XT::Grid::is_layer<GridLayerType>::value,
    std::unique_ptr<
        WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight,
                                 MatrixType& matrix,
                                 const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>(
      over_integrate, weight, matrix, range_space, source_space, grid_layer);
}


// ////////////////// //
// WeightedL2Operator //
// ////////////////// //

// forward, needed for the traits
template <class WeightFunctionType, class GridLayer, class Field = typename WeightFunctionType::RangeFieldType>
class WeightedL2Operator;


namespace internal {


template <class WeightFunctionType, class GridLayerType, class Field>
class WeightedL2OperatorTraits
{
public:
  typedef WeightedL2Operator<WeightFunctionType, GridLayerType, Field> derived_type;
  typedef NoJacobian JacobianType;
  typedef Field FieldType;
};


} // namespace internal


template <class WeightFunctionType, class GridLayerType, class Field>
class WeightedL2Operator
  : public OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridLayerType, Field>>
{
  typedef OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridLayerType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  WeightedL2Operator(const WeightFunctionType& weight, GridLayerType grid_layer, const size_t over_integrate = 0)
    : weight_(weight)
    , grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {}

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range,
             const XT::Common::Parameter& param = {}) const
  {
    typedef
        typename XT::LA::Container<typename VectorType::ScalarType, VectorType::Traits::sparse_matrix_type>::MatrixType
            MatrixType;
    auto op = make_weighted_l2_matrix_operator<MatrixType>(
        weight_, source.space(), range.space(), grid_layer_, over_integrate_);
    op->apply(source, range, param);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
                   const XT::Common::Parameter& param = {}) const
  {
    auto product = make_weighted_l2_localizable_product(weight_, grid_layer_, range, source, over_integrate_, param);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void
  apply_inverse(const RangeType& /*range*/, SourceType& /*source*/, const XT::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "yet");
    return {"depends_on_the_vector_type_of_the_discrete_function"};
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

private:
  const WeightFunctionType& weight_;
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class WeightedL2Operator


// ///////////////////////// //
// make_weighted_l2_operator //
// ///////////////////////// //

template <class GridLayerType, class WeightFunctionType>
typename std::enable_if<
    XT::Functions::is_localizable_function<WeightFunctionType>::value && XT::Grid::is_layer<GridLayerType>::value,
    std::unique_ptr<
        WeightedL2Operator<WeightFunctionType, GridLayerType, typename WeightFunctionType::RangeFieldType>>>::type
make_weighted_l2_operator(const GridLayerType& grid_layer,
                          const WeightFunctionType& weight,
                          const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      WeightedL2Operator<WeightFunctionType, GridLayerType, typename WeightFunctionType::RangeFieldType>>(
      weight, grid_layer, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
