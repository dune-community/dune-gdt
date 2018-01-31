// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_LAPLACE_HH
#define DUNE_GDT_OPERATORS_LAPLACE_HH

#include <limits>
#include <type_traits>

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/type_traits.hh>

#include "elliptic.hh"

namespace Dune {
namespace GDT {


// ///////////////////////// //
// LaplaceLocalizableProduct //
// ///////////////////////// //

template <class GridLayer, class Range, class Source = Range, class Field = typename Range::RangeFieldType>
class LaplaceLocalizableProduct
    : XT::Common::ConstStorageProvider<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                       typename GridLayer::ctype,
                                                                       GridLayer::dimension,
                                                                       Field,
                                                                       1>>,
      public EllipticLocalizableProduct<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                        typename GridLayer::ctype,
                                                                        GridLayer::dimension,
                                                                        Field,
                                                                        1>,
                                        void,
                                        GridLayer,
                                        Range,
                                        Source,
                                        Field>
{
  typedef XT::Common::ConstStorageProvider<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                           typename GridLayer::ctype,
                                                                           GridLayer::dimension,
                                                                           Field,
                                                                           1>>
      FunctionProvider;
  typedef EllipticLocalizableProduct<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                     typename GridLayer::ctype,
                                                                     GridLayer::dimension,
                                                                     Field,
                                                                     1>,
                                     void,
                                     GridLayer,
                                     Range,
                                     Source,
                                     Field>
      BaseType;

  // We suffer from the same problem as in L2LocalizableProduct, see there for an explanation.
  template <bool anything>
  struct tag
  {
    explicit tag(int)
    {
    }
  };

  template <class... Args>
  explicit LaplaceLocalizableProduct(tag<false>, Args&&... args)
    : FunctionProvider(FunctionProvider::make(1.))
    , BaseType(FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

  template <class... Args>
  explicit LaplaceLocalizableProduct(tag<true>, const size_t over_integrate, Args&&... args)
    : FunctionProvider(FunctionProvider::make(1.))
    , BaseType(over_integrate, FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

public:
  /**
   * \brief Creates a localizable Laplace product (aka a localizable semi-H1 product).
   *
   *        We suffer from the same problems as L2LocalizableProduct, see also the documentation of
   *        \sa L2LocalizableProduct::L2LocalizableProduct(). This ctor can be used as follows, where over_integrate is
   *        a non-negative integer and ...args are the arguments for LocalizableProductBase (i.e., a grid_layer, a
   *        range and possibly a source):
\code
LaplaceLocalizableProduct(over_integrate, ...args);
LaplaceLocalizableProduct(...args);
\endcode
   */
  template <typename possibly_int_t,
            class... Args,
            typename std::enable_if<!std::is_same<possibly_int_t, tag<true>>::value
                                        && !std::is_same<possibly_int_t, tag<false>>::value,
                                    int>::type = 0>
  explicit LaplaceLocalizableProduct(possibly_int_t&& possibly_over_integrate, Args&&... args)
    : LaplaceLocalizableProduct(tag<std::numeric_limits<typename std::decay<possibly_int_t>::type>::is_integer>(0),
                                std::forward<possibly_int_t>(possibly_over_integrate),
                                std::forward<Args>(args)...)
  {
  }
}; // class LaplaceLocalizableProduct


// //////////////////////////////// //
// make_laplace_localizable_product //
// //////////////////////////////// //

/**
 * \sa LaplaceLocalizableProduct
 */
template <class GridLayerType, class RangeType, class SourceType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<RangeType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value,
                        std::unique_ptr<LaplaceLocalizableProduct<GridLayerType, RangeType, SourceType>>>::type
make_laplace_localizable_product(const GridLayerType& grid_layer,
                                 const RangeType& range,
                                 const SourceType& source,
                                 const size_t over_integrate = 0,
                                 const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<LaplaceLocalizableProduct<GridLayerType, RangeType, SourceType>>(
      over_integrate, grid_layer, range, source, param);
}


// ///////////////////// //
// LaplaceMatrixOperator //
// ///////////////////// //

template <class RangeSpace,
          class Matrix = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridLayer = typename RangeSpace::GridLayerType,
          class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class LaplaceMatrixOperator
    : XT::Common::ConstStorageProvider<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                       typename GridLayer::ctype,
                                                                       GridLayer::dimension,
                                                                       Field,
                                                                       1>>,
      public EllipticMatrixOperator<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                    typename GridLayer::ctype,
                                                                    GridLayer::dimension,
                                                                    Field,
                                                                    1>,
                                    void,
                                    RangeSpace,
                                    Matrix,
                                    GridLayer,
                                    SourceSpace,
                                    Field>
{
  typedef XT::Common::ConstStorageProvider<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                           typename GridLayer::ctype,
                                                                           GridLayer::dimension,
                                                                           Field,
                                                                           1>>
      FunctionProvider;
  typedef EllipticMatrixOperator<XT::Functions::ConstantFunction<XT::Grid::extract_entity_t<GridLayer>,
                                                                 typename GridLayer::ctype,
                                                                 GridLayer::dimension,
                                                                 Field,
                                                                 1>,
                                 void,
                                 RangeSpace,
                                 Matrix,
                                 GridLayer,
                                 SourceSpace,
                                 Field>
      BaseType;

  // We suffer from the same problem as in L2LocalizableProduct, see there for an explanation.
  template <bool anything>
  struct tag
  {
    explicit tag(int)
    {
    }
  };

  template <class... Args>
  explicit LaplaceMatrixOperator(tag<false>, Args&&... args)
    : FunctionProvider(FunctionProvider::make(1.))
    , BaseType(FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

  template <class... Args>
  explicit LaplaceMatrixOperator(tag<true>, const size_t over_integrate, Args&&... args)
    : FunctionProvider(FunctionProvider::make(1.))
    , BaseType(over_integrate, FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

public:
  /**
   * \brief Creates a matrix-based Laplace operator.
   *
   *        We suffer from the same problems as L2LocalizableProduct, see also the documentation of
   *        \sa L2LocalizableProduct::L2LocalizableProduct(). This ctor can be used as follows, where over_integrate is
   *        a non-negative integer and ...args are the arguments for MatrixOperatorBase (i.e., possibly a matrix,
   *        a range space, possibly a grid_layer and possibly a source space):
\code
LaplaceLocalizableProduct(over_integrate, ...args);
LaplaceLocalizableProduct(...args);
\endcode
   *        If no matrix is provided, an appropriate matrix of given MatrixType will be created and is accessible via
   *        matrix().
   */
  template <typename possibly_int_t,
            class... Args,
            typename std::enable_if<!std::is_same<possibly_int_t, tag<true>>::value
                                        && !std::is_same<possibly_int_t, tag<false>>::value,
                                    int>::type = 0>
  explicit LaplaceMatrixOperator(possibly_int_t&& possibly_over_integrate, Args&&... args)
    : LaplaceMatrixOperator(tag<std::numeric_limits<typename std::decay<possibly_int_t>::type>::is_integer>(0),
                            std::forward<possibly_int_t>(possibly_over_integrate),
                            std::forward<Args>(args)...)
  {
  }
}; // class LaplaceMatrixOperator


// /////////////////////// //
// make_laplace_matrix_operator //
// /////////////////////// //

// without matrix

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(space);
\endcode
 */
template <class MatrixType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType>>>::type
make_laplace_matrix_operator(const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType>>(over_integrate, space);
}

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(space, grid_layer);
\endcode
 */
template <class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType, GridLayerType>>>::type
make_laplace_matrix_operator(const SpaceType& space, const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType, GridLayerType>>(
      over_integrate, space, grid_layer);
}

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(range_space, source_space, grid_layer);
\endcode
 */
template <class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>>::type
    make_laplace_matrix_operator(const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>(
          over_integrate, range_space, source_space, grid_layer);
}

// with matrix

/**
 * \brief Creates a Laplace matrix operator (source and range space are given by space, grid_layer of the space is
 * used).
 */
template <class MatrixType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType>>>::type
make_laplace_matrix_operator(MatrixType& matrix, const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType>>(over_integrate, matrix, space);
}

/**
 * \brief Creates a Laplace matrix operator (source and range space are given by space).
 */
template <class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType, GridLayerType>>>::type
make_laplace_matrix_operator(MatrixType& matrix,
                             const SpaceType& space,
                             const GridLayerType& grid_layer,
                             const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType, GridLayerType>>(
      over_integrate, matrix, space, grid_layer);
}

/**
 * \brief Creates a Laplace matrix operator.
 */
template <class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>>::type
    make_laplace_matrix_operator(MatrixType& matrix,
                                 const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space,
                                 const GridLayerType& grid_layer,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridLayerType, SourceSpaceType>>(
          over_integrate, matrix, range_space, source_space, grid_layer);
}


// /////////////// //
// LaplaceOperator //
// /////////////// //

// forward, needed for the traits
template <class GridLayer, class Field = double>
class LaplaceOperator;


namespace internal {


template <class GridLayerType, class Field>
class LaplaceOperatorTraits
{
public:
  typedef LaplaceOperator<GridLayerType, Field> derived_type;
  typedef NoJacobian JacobianType;
  typedef Field FieldType;
};


} // namespace internal


template <class GridLayerType, class Field>
class LaplaceOperator : public OperatorInterface<internal::LaplaceOperatorTraits<GridLayerType, Field>>
{
  typedef OperatorInterface<internal::LaplaceOperatorTraits<GridLayerType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  template <class... Args>
  LaplaceOperator(GridLayerType grid_layer, const size_t over_integrate = 0)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range,
             const XT::Common::Parameter& param = {}) const
  {
    typedef
        typename XT::LA::Container<typename VectorType::ScalarType, VectorType::Traits::sparse_matrix_type>::MatrixType
            MatrixType;
    auto op = make_laplace_matrix_operator<MatrixType>(source.space(), range.space(), grid_layer_, over_integrate_);
    op->apply(source, range, param);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
                   const XT::Common::Parameter& param = {}) const
  {
    auto product = make_laplace_localizable_product(grid_layer_, range, source, over_integrate_, param);
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
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class LaplaceOperator


// ///////////////////// //
// make_laplace_operator //
// ///////////////////// //

template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value, std::unique_ptr<LaplaceOperator<GridLayerType>>>::type
make_laplace_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<LaplaceOperator<GridLayerType>>(grid_layer, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LAPLACE_HH
