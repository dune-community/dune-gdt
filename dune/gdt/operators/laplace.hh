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

template< class GridView,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType >
class LaplaceLocalizableProduct
  : Stuff::Common::ConstStorageProvider<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type,
                                    typename GridView::ctype, GridView::dimension, Field, 1 > >
  , public EllipticLocalizableProduct<
        Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 >,
        void, GridView, Range, Source, Field >
{
  typedef Stuff::Common::ConstStorageProvider
      < Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                    GridView::dimension, Field, 1 > > FunctionProvider;
  typedef EllipticLocalizableProduct<
      Stuff::Functions::Constant< typename Stuff::Grid::Entity< GridView >::type, typename GridView::ctype,
                                  GridView::dimension, Field, 1 >,
      void, GridView, Range, Source, Field > BaseType;

  // We suffer from the same problem as in L2LocalizableProduct, see there for an explanation.
  template< bool anything >
  struct tag{ explicit tag(int) {}};

  template< class ...Args >
  explicit LaplaceLocalizableProduct(tag< false >, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(FunctionProvider::access(), std::forward< Args >(args)...)
  {}

  template< class ...Args >
  explicit LaplaceLocalizableProduct(tag< true >, const size_t over_integrate, Args&& ...args)
    : FunctionProvider(1.)
    , BaseType(over_integrate, FunctionProvider::access(), std::forward< Args >(args)...)
  {}

public:
  /**
   * \brief Creates a localizable Laplace product (aka a localizable semi-H1 product).
   *
   *        We suffer from the same problems as L2LocalizableProduct, see also the documentation of
   *        \sa L2LocalizableProduct::L2LocalizableProduct(). This ctor can be used as follows, where over_integrate is
   *        a non-negative integer and ...args are the arguments for LocalizableProductDefault (i.e., a grid_view, a
   *        range and possibly a source):
\code
LaplaceLocalizableProduct(over_integrate, ...args);
LaplaceLocalizableProduct(...args);
\endcode
   */
  template< typename possibly_int_t,
            class ...Args,
            typename std::enable_if< !std::is_same< possibly_int_t, tag< true > >::value, int >::type = 0 >
  explicit LaplaceLocalizableProduct(possibly_int_t&& possibly_over_integrate, Args&& ...args)
    : LaplaceLocalizableProduct(tag< std::numeric_limits< typename std::decay< possibly_int_t >::type >::is_integer >(0),
                                std::forward< possibly_int_t >(possibly_over_integrate),
                                std::forward< Args >(args)...)
  {}
}; // class LaplaceLocalizableProduct


// //////////////////////////////// //
// make_laplace_localizable_product //
// //////////////////////////////// //

/**
 * \sa LaplaceLocalizableProduct
 */
template< class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                            , std::unique_ptr< LaplaceLocalizableProduct< GridViewType, RangeType, SourceType > >
                            >::type
make_laplace_localizable_product(const GridViewType& grid_view,
                                 const RangeType& range,
                                 const SourceType& source,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique< LaplaceLocalizableProduct< GridViewType, RangeType, SourceType > >(
        over_integrate, grid_view, range, source);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LAPLACE_HH
