// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_OPERATORS_WEIGHTED_L2_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


// //////////////////////////// //
// WeightedL2LocalizableProduct //
// //////////////////////////// //

template <class WeightFunctionType, class GridView, class Range, class Source = Range,
          class Field = typename Range::RangeFieldType>
class WeightedL2LocalizableProduct : public LocalizableProductBase<GridView, Range, Source, Field>
{
  typedef LocalizableProductBase<GridView, Range, Source, Field> BaseType;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>> LocalWeightedL2OperatorType;

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
                            && XT::Grid::is_layer<GridViewType>::value
                            && Stuff::is_localizable_function<RangeType>::value
                            && Stuff::is_localizable_function<SourceType>::value,
                        std::unique_ptr<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType,
                                                                     SourceType>>>::type
make_weighted_l2_localizable_product(const WeightFunctionType& weight, const GridViewType& grid_view,
                                     const RangeType& range, const SourceType& source, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType, SourceType>>(
          over_integrate, weight, grid_view, range, source);
}


// //////////////////////// //
// WeightedL2MatrixOperator //
// //////////////////////// //

template <class WeightFunctionType, class RangeSpace,
          class Matrix   = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridView = typename RangeSpace::GridViewType, class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class WeightedL2MatrixOperator
    : public MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume>
{
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume> BaseType;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>> LocalWeightedL2OperatorType;

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
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
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
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space, grid_view);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType, class GridViewType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && Stuff::is_localizable_function<WeightFunctionType>::value
                  && is_space<SpaceType>::value && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
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
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value
                            && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space, const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
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
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, matrix, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType, class GridViewType>
typename std::
    enable_if<Stuff::is_localizable_function<WeightFunctionType>::value && XT::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
          over_integrate, weight, matrix, space, grid_view);
}

/**
 * \brief Creates a weighted L2 matrix operator.
 */
template <class WeightFunctionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix,
                                 const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                 const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
                                                                RangeSpaceType,
                                                                MatrixType,
                                                                GridViewType,
                                                                SourceSpaceType>>(
      over_integrate, weight, matrix, range_space, source_space, grid_view);
}


// ////////////////// //
// WeightedL2Operator //
// ////////////////// //

// forward, needed for the traits
template <class WeightFunctionType, class GridView, class Field = typename WeightFunctionType::RangeFieldType>
class WeightedL2Operator;


namespace internal {


template <class WeightFunctionType, class GridViewType, class Field>
class WeightedL2OperatorTraits
{
public:
  typedef WeightedL2Operator<WeightFunctionType, GridViewType, Field> derived_type;
  typedef Field FieldType;
};


} // namespace internal


template <class WeightFunctionType, class GridViewType, class Field>
class WeightedL2Operator
    : public OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridViewType, Field>>
{
  typedef OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridViewType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  WeightedL2Operator(const WeightFunctionType& weight, GridViewType grid_view, const size_t over_integrate = 0)
    : weight_(weight)
    , grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range) const
  {
    typedef typename XT::LA::Container<typename VectorType::ScalarType, VectorType::Traits::sparse_matrix_type>::MatrixType
        MatrixType;
    auto op = make_weighted_l2_matrix_operator<MatrixType>(
        weight_, source.space(), range.space(), grid_view_, over_integrate_);
    op->apply(source, range);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& source) const
  {
    auto product = make_weighted_l2_localizable_product(weight_, grid_view_, range, source, over_integrate_);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/, SourceType& /*source*/,
                     const XT::Common::Configuration& /*opts*/) const
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
  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class WeightedL2Operator


// ///////////////////////// //
// make_weighted_l2_operator //
// ///////////////////////// //

template <class GridViewType, class WeightFunctionType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2Operator<WeightFunctionType, GridViewType,
                                                           typename WeightFunctionType::RangeFieldType>>>::type
make_weighted_l2_operator(const GridViewType& grid_view, const WeightFunctionType& weight,
                          const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<WeightedL2Operator<WeightFunctionType, GridViewType, typename WeightFunctionType::RangeFieldType>>(
          weight, grid_view, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
