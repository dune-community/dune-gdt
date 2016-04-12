// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_LAPLACE_HH
#define DUNE_GDT_OPERATORS_LAPLACE_HH

#include <limits>
#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/grid/entity.hh>

#include "elliptic.hh"

namespace Dune {
namespace GDT {


// ///////////////////////// //
// LaplaceLocalizableProduct //
// ///////////////////////// //

template <class GridView, class Range, class Source = Range, class Field = typename Range::RangeFieldType>
class LaplaceLocalizableProduct
    : Stuff::Common::ConstStorageProvider<Stuff::Functions::Constant<
          typename Stuff::Grid::Entity<GridView>::type, typename GridView::ctype, GridView::dimension, Field, 1>>,
      public EllipticLocalizableProduct<Stuff::Functions::Constant<typename Stuff::Grid::Entity<GridView>::type,
                                                                   typename GridView::ctype, GridView::dimension, Field,
                                                                   1>,
                                        void, GridView, Range, Source, Field>
{
  typedef Stuff::Common::ConstStorageProvider<Stuff::Functions::Constant<
      typename Stuff::Grid::Entity<GridView>::type, typename GridView::ctype, GridView::dimension, Field, 1>>
      FunctionProvider;
  typedef EllipticLocalizableProduct<Stuff::Functions::Constant<typename Stuff::Grid::Entity<GridView>::type,
                                                                typename GridView::ctype, GridView::dimension, Field,
                                                                1>,
                                     void, GridView, Range, Source, Field> BaseType;

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
    : FunctionProvider(1.)
    , BaseType(FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

  template <class... Args>
  explicit LaplaceLocalizableProduct(tag<true>, const size_t over_integrate, Args&&... args)
    : FunctionProvider(1.)
    , BaseType(over_integrate, FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

public:
  /**
   * \brief Creates a localizable Laplace product (aka a localizable semi-H1 product).
   *
   *        We suffer from the same problems as L2LocalizableProduct, see also the documentation of
   *        \sa L2LocalizableProduct::L2LocalizableProduct(). This ctor can be used as follows, where over_integrate is
   *        a non-negative integer and ...args are the arguments for LocalizableProductBase (i.e., a grid_view, a
   *        range and possibly a source):
\code
LaplaceLocalizableProduct(over_integrate, ...args);
LaplaceLocalizableProduct(...args);
\endcode
   */
  template <typename possibly_int_t, class... Args,
            typename std::enable_if<!std::is_same<possibly_int_t, tag<true>>::value
                                        && !std::is_same<possibly_int_t, tag<false>>::value,
                                    int>::type = 0>
  explicit LaplaceLocalizableProduct(possibly_int_t&& possibly_over_integrate, Args&&... args)
    : LaplaceLocalizableProduct(tag<std::numeric_limits<typename std::decay<possibly_int_t>::type>::is_integer>(0),
                                std::forward<possibly_int_t>(possibly_over_integrate), std::forward<Args>(args)...)
  {
  }
}; // class LaplaceLocalizableProduct


// //////////////////////////////// //
// make_laplace_localizable_product //
// //////////////////////////////// //

/**
 * \sa LaplaceLocalizableProduct
 */
template <class GridViewType, class RangeType, class SourceType>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value
                            && Stuff::is_localizable_function<RangeType>::value
                            && Stuff::is_localizable_function<SourceType>::value,
                        std::unique_ptr<LaplaceLocalizableProduct<GridViewType, RangeType, SourceType>>>::type
make_laplace_localizable_product(const GridViewType& grid_view, const RangeType& range, const SourceType& source,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceLocalizableProduct<GridViewType, RangeType, SourceType>>(
      over_integrate, grid_view, range, source);
}


// ///////////////////// //
// LaplaceMatrixOperator //
// ///////////////////// //

template <class RangeSpace,
          class Matrix   = typename Stuff::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridView = typename RangeSpace::GridViewType, class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class LaplaceMatrixOperator
    : Stuff::Common::ConstStorageProvider<Stuff::Functions::Constant<
          typename Stuff::Grid::Entity<GridView>::type, typename GridView::ctype, GridView::dimension, Field, 1>>,
      public EllipticMatrixOperator<Stuff::Functions::Constant<typename Stuff::Grid::Entity<GridView>::type,
                                                               typename GridView::ctype, GridView::dimension, Field, 1>,
                                    void, RangeSpace, Matrix, GridView, SourceSpace, Field>
{
  typedef Stuff::Common::ConstStorageProvider<Stuff::Functions::Constant<
      typename Stuff::Grid::Entity<GridView>::type, typename GridView::ctype, GridView::dimension, Field, 1>>
      FunctionProvider;
  typedef EllipticMatrixOperator<Stuff::Functions::Constant<typename Stuff::Grid::Entity<GridView>::type,
                                                            typename GridView::ctype, GridView::dimension, Field, 1>,
                                 void, RangeSpace, Matrix, GridView, SourceSpace, Field> BaseType;

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
    : FunctionProvider(1.)
    , BaseType(FunctionProvider::access(), std::forward<Args>(args)...)
  {
  }

  template <class... Args>
  explicit LaplaceMatrixOperator(tag<true>, const size_t over_integrate, Args&&... args)
    : FunctionProvider(1.)
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
   *        a range space, possibly a grid_view and possibly a source space):
\code
LaplaceLocalizableProduct(over_integrate, ...args);
LaplaceLocalizableProduct(...args);
\endcode
   *        If no matrix is provided, an appropriate matrix of given MatrixType will be created and is accessible via
   *        matrix().
   */
  template <typename possibly_int_t, class... Args,
            typename std::enable_if<!std::is_same<possibly_int_t, tag<true>>::value
                                        && !std::is_same<possibly_int_t, tag<false>>::value,
                                    int>::type = 0>
  explicit LaplaceMatrixOperator(possibly_int_t&& possibly_over_integrate, Args&&... args)
    : LaplaceMatrixOperator(tag<std::numeric_limits<typename std::decay<possibly_int_t>::type>::is_integer>(0),
                            std::forward<possibly_int_t>(possibly_over_integrate), std::forward<Args>(args)...)
  {
  }
}; // class LaplaceMatrixOperator


// /////////////////////// //
// make_laplace_matrix_operator //
// /////////////////////// //

// without matrix

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(space);
\endcode
 */
template <class MatrixType, class SpaceType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType>>>::type
make_laplace_matrix_operator(const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType>>(over_integrate, space);
}

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(space, grid_view);
\endcode
 */
template <class MatrixType, class SpaceType, class GridViewType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType>>>::type
make_laplace_matrix_operator(const SpaceType& space, const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType>>(over_integrate, space, grid_view);
}

/**
 * \brief Creates a Laplace matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_laplace_matrix_operator< MatrixType >(range_space, source_space, grid_view);
\endcode
 */
template <class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::
    enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridViewType, SourceSpaceType>>>::type
    make_laplace_matrix_operator(const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                 const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridViewType, SourceSpaceType>>(
      over_integrate, range_space, source_space, grid_view);
}

// with matrix

/**
 * \brief Creates a Laplace matrix operator (source and range space are given by space, grid_view of the space is used).
 */
template <class MatrixType, class SpaceType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType>>>::type
make_laplace_matrix_operator(MatrixType& matrix, const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType>>(over_integrate, matrix, space);
}

/**
 * \brief Creates a Laplace matrix operator (source and range space are given by space).
 */
template <class MatrixType, class SpaceType, class GridViewType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType>>>::type
make_laplace_matrix_operator(MatrixType& matrix, const SpaceType& space, const GridViewType& grid_view,
                             const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<SpaceType, MatrixType, GridViewType>>(
      over_integrate, matrix, space, grid_view);
}

/**
 * \brief Creates a Laplace matrix operator.
 */
template <class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::
    enable_if<Stuff::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridViewType, SourceSpaceType>>>::type
    make_laplace_matrix_operator(MatrixType& matrix, const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space, const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceMatrixOperator<RangeSpaceType, MatrixType, GridViewType, SourceSpaceType>>(
      over_integrate, matrix, range_space, source_space, grid_view);
}


// /////////////// //
// LaplaceOperator //
// /////////////// //

// forward, needed for the traits
template <class GridView, class Field = double>
class LaplaceOperator;


namespace internal {


template <class GridViewType, class Field>
class LaplaceOperatorTraits
{
public:
  typedef LaplaceOperator<GridViewType, Field> derived_type;
  typedef Field FieldType;
};


} // namespace internal


template <class GridViewType, class Field>
class LaplaceOperator : public OperatorInterface<internal::LaplaceOperatorTraits<GridViewType, Field>>
{
  typedef OperatorInterface<internal::LaplaceOperatorTraits<GridViewType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  template <class... Args>
  LaplaceOperator(GridViewType grid_view, const size_t over_integrate = 0)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range) const
  {
    typedef typename Stuff::LA::Container<typename VectorType::ScalarType, VectorType::sparse_matrix_type>::MatrixType
        MatrixType;
    auto op = make_laplace_matrix_operator<MatrixType>(source.space(), range.space(), grid_view_, over_integrate_);
    op->apply(source, range);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& source) const
  {
    auto product = make_laplace_localizable_product(grid_view_, range, source, over_integrate_);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/, SourceType& /*source*/,
                     const Stuff::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "yet");
    return {"depends_on_the_vector_type_of_the_discrete_function"};
  }

  Stuff::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

private:
  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class LaplaceOperator


// ///////////////////// //
// make_laplace_operator //
// ///////////////////// //

template <class GridViewType>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<LaplaceOperator<GridViewType>>>::type
make_laplace_operator(const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<LaplaceOperator<GridViewType>>(grid_view, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LAPLACE_HH
