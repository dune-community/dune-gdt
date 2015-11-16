// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_L2_HH
#define DUNE_GDT_OPERATORS_L2_HH

#include <limits>
#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/grid/entity.hh>

#include "weighted-l2.hh"

namespace Dune {
namespace GDT {


// //////////////////// //
// L2LocalizableProduct //
// //////////////////// //

template< class GridView,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType >
class L2LocalizableProduct
  : Stuff::Common::ConstStorageProvider<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type,
                                    typename GridView::ctype, GridView::dimension, Field, 1 > >
  , public WeightedL2LocalizableProduct<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 >,
        GridView, Range, Source, Field >
{
  typedef Stuff::Common::ConstStorageProvider
      < Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 > > FunctionProvider;
  typedef WeightedL2LocalizableProduct<
      Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                  GridView::dimension, Field, 1 >,
      GridView, Range, Source, Field > BaseType;

  // The following tag and the two ctors are unfortunately necessary. There should have been two ctors,
  //     L2LocalizableProduct(Args&& ...args)
  // and
  //     L2LocalizableProduct(const size_t over_integrate, Args&& ...args),
  // both using perfect forwarding on ...args (see the documentation of the ctor below for application szenarios).
  // However, the compiler was not able to differentiate between the two. The work-around below was inspired by:
  //     http://stackoverflow.com/questions/32957830/c-variadic-template-constructor-and-common-constructors
  template< bool anything >
  struct tag{ explicit tag(int) {}}; // we do not want tag to be implicitely constructable, just to be safe

  template< class ...Args >
  explicit L2LocalizableProduct(tag< false >, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(FunctionProvider::access(), std::forward< Args >(args)...)
  {}

  template< class ...Args >
  explicit L2LocalizableProduct(tag< true >, const size_t over_integrate, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(over_integrate, FunctionProvider::access(), std::forward< Args >(args)...)
  {}

public:
  /**
   * \brief Creates a localizable L2 product.
   *
   *        Since the L2LocalizableProduct does not require anything beyond what LocalizableProductDefault does,
   *        we forward all arguments to WeightedL2LocalizableProduct (which, in turn, forwards them to
   *        LocalizableProductDefault). In addition, we want to allow to specify an over_integrate parameter, which is
   *        passen on to WeightedL2LocalizableProduct. This leads to two kind of ctors, where over_integrate is a
   *        non-negative integer and ...args are the arguments for LocalizableProductDefault (i.e., a grid_view, a range
   *        and possibly a source):
\code
L2LocalizableProduct(over_integrate, ...args);
L2LocalizableProduct(...args);
\endcode
   *        For technical reasons we require a work-around to realize these ctors, leading to the present weird
   *        signature. Nevertheless, you can just use it as in the example above.
   */
  template< typename possibly_int_t,
            class ...Args,
            typename std::enable_if< !std::is_same< possibly_int_t, tag< true > >::value, int >::type = 0 >
  explicit L2LocalizableProduct(possibly_int_t&& possibly_over_integrate, Args&& ...args)
    : L2LocalizableProduct(tag< std::numeric_limits< typename std::decay< possibly_int_t >::type >::is_integer >(0),
                           std::forward< possibly_int_t >(possibly_over_integrate),
                           std::forward< Args >(args)...)
  {}
}; // class L2LocalizableProduct


// /////////////////////////// //
// make_l2_localizable_product //
// /////////////////////////// //

/**
 * \sa L2LocalizableProduct
 */
template< class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                            , std::unique_ptr< L2LocalizableProduct< GridViewType, RangeType, SourceType > >
                            >::type
make_l2_localizable_product(const GridViewType& grid_view,
                            const RangeType& range,
                            const SourceType& source,
                            const size_t over_integrate = 0)
{
  return DSC::make_unique< L2LocalizableProduct< GridViewType, RangeType, SourceType > >(
        over_integrate, grid_view, range, source);
}


// //////////////// //
// L2MatrixOperator //
// //////////////// //

template< class RangeSpace,
          class Matrix = typename Stuff::LA::Container< typename RangeSpace::RangeFieldType >::MatrixType,
          class GridView = typename RangeSpace::GridViewType,
          class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType >
class L2MatrixOperator
  : Stuff::Common::ConstStorageProvider<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type,
                                    typename GridView::ctype, GridView::dimension, Field, 1 > >
  , public WeightedL2MatrixOperator<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 >,
        RangeSpace, Matrix, GridView, SourceSpace, Field >
{
  typedef Stuff::Common::ConstStorageProvider
      < Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 > > FunctionProvider;
  typedef WeightedL2MatrixOperator<
      Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                  GridView::dimension, Field, 1 >,
      RangeSpace, Matrix, GridView, SourceSpace, Field >BaseType;

  // We suffer from the same problem as in L2LocalizableProduct, see above for an explanation.
  template< bool anything >
  struct tag{ explicit tag(int) {}};

  template< class ...Args >
  explicit L2MatrixOperator(tag< false >, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(FunctionProvider::access(), std::forward< Args >(args)...)
  {}

  template< class ...Args >
  explicit L2MatrixOperator(tag< true >, const size_t over_integrate, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(over_integrate, FunctionProvider::access(), std::forward< Args >(args)...)
  {}

public:
  /**
   * \brief Creates a matrix-based L2 operator.
   *
   *        We suffer from the same problems as L2LocalizableProduct, see also the documentation of
   *        \sa L2LocalizableProduct::L2LocalizableProduct(). This ctor can be used as follows, where over_integrate is
   *        a non-negative integer and ...args are the arguments for MatrixOperatorDefault (i.e., possibly a matrix,
   *        a range space, possibly a grid_view and possibly a source space):
\code
L2LocalizableProduct(over_integrate, ...args);
L2LocalizableProduct(...args);
\endcode
   *        If no matrix is provided, an appropriate matrix of given MatrixType will be created and is accessible via
   *        matrix().
   */
  template< typename possibly_int_t,
            class ...Args,
            typename std::enable_if< !std::is_same< possibly_int_t, tag< true > >::value, int >::type = 0 >
  explicit L2MatrixOperator(possibly_int_t&& possibly_over_integrate, Args&& ...args)
    : L2MatrixOperator(tag< std::numeric_limits< typename std::decay< possibly_int_t >::type >::is_integer >(0),
                       std::forward< possibly_int_t >(possibly_over_integrate),
                       std::forward< Args >(args)...)
  {}
}; // class L2MatrixOperator


// /////////////////////// //
// make_l2_matrix_operator //
// /////////////////////// //

// without matrix

/**
 * \brief Creates an L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically, source and
 *        range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_l2_matrix_operator< MatrixType >(space);
\endcode
 */
template< class MatrixType, class SpaceType >
    typename std::enable_if< Stuff::LA::is_matrix< MatrixType >::value && is_space< SpaceType >::value
                           , std::unique_ptr< L2MatrixOperator< SpaceType, MatrixType > >
                           >::type
make_l2_matrix_operator(const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< SpaceType, MatrixType > >(over_integrate, space);
}

/**
 * \brief Creates an L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically, source and
 *        range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_l2_matrix_operator< MatrixType >(space, grid_view);
\endcode
 */
template< class MatrixType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2MatrixOperator< SpaceType, MatrixType, GridViewType > >
                           >::type
make_l2_matrix_operator(const SpaceType& space, const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< SpaceType, MatrixType, GridViewType > >(over_integrate, space, grid_view);
}

/**
 * \brief Creates an L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_l2_matrix_operator< MatrixType >(range_space, source_space, grid_view);
\endcode
 */
template< class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2MatrixOperator< RangeSpaceType, MatrixType, GridViewType,
                                                                SourceSpaceType > >
                           >::type
make_l2_matrix_operator(const RangeSpaceType& range_space,
                        const SourceSpaceType& source_space,
                        const GridViewType& grid_view,
                        const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
      (over_integrate, range_space, source_space, grid_view);
}

// with matrix

/**
 * \brief Creates an L2 matrix operator (source and range space are given by space, grid_view of the space is
 *        used).
 */
template< class MatrixType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2MatrixOperator< SpaceType, MatrixType > >
                           >::type
make_l2_matrix_operator(MatrixType& matrix, const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< SpaceType, MatrixType > >(over_integrate, matrix, space);
}

/**
 * \brief Creates an L2 matrix operator (source and range space are given by space).
 */
template< class MatrixType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2MatrixOperator< SpaceType, MatrixType, GridViewType > >
                           >::type
make_l2_matrix_operator(MatrixType& matrix,
                        const SpaceType& space,
                        const GridViewType& grid_view,
                        const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< SpaceType, MatrixType, GridViewType > >
      (over_integrate, matrix, space, grid_view);
}

/**
 * \brief Creates an L2 matrix operator.
 */
template< class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2MatrixOperator< RangeSpaceType, MatrixType, GridViewType,
                                                                SourceSpaceType > >
                           >::type
make_l2_matrix_operator(MatrixType& matrix,
                        const RangeSpaceType& range_space,
                        const SourceSpaceType& source_space,
                        const GridViewType& grid_view,
                        const size_t over_integrate = 0)
{
  return DSC::make_unique< L2MatrixOperator< RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
      (over_integrate, matrix, range_space, source_space, grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_L2_HH
