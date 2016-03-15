// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_HH

#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/integrals.hh>
#include <dune/gdt/spaces/interface.hh>

#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace GDT {


// ////////////////////////// //
// EllipticLocalizableProduct //
// ////////////////////////// //

template< class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class GridView,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType >
class EllipticLocalizableProduct
  : public LocalizableProductDefault< GridView, Range, Source, Field >
{
  typedef LocalizableProductDefault< GridView, Range, Source, Field >           BaseType;
  typedef LocalVolumeIntegralOperator
      < LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > > LocalEllipticOperatorType;
public:
  // Usually, we only have to hold the data functions for the local operator and perfect forward the rest of the
  // arguments to BaseType. Here it is a bit more complicated, since DiffusionTensorType might be void (and
  // DiffusionFactorType is the the only diffusion). To handle this case we require the enable_if hacks below to
  // disable half of the ctors if DiffusionTensorType is void. In addition we require each ctor twice, once with
  // over_integrate and once without (since the perfect forwarding does not allow for default arguments).

  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value) //
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const DiffusionImp& diffusion, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const size_t over_integrate,
                                      const DiffusionImp& diffusion,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const size_t over_integrate,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

private:
  const LocalEllipticOperatorType local_elliptic_operator_;
}; // class EllipticLocalizableProduct


// ///////////////////////////////// //
// make_elliptic_localizable_product //
// ///////////////////////////////// //

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion.
 */
template< class DiffusionType, class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                            , std::unique_ptr< EllipticLocalizableProduct< DiffusionType, void, GridViewType,
                                                                           RangeType, SourceType > >
                            >::type
make_elliptic_localizable_product(const DiffusionType& diffusion,
                                  const GridViewType& grid_view,
                                  const RangeType& range,
                                  const SourceType& source,
                                  const size_t over_integrate = 0)
{
  return DSC::make_unique
      < EllipticLocalizableProduct< DiffusionType, void, GridViewType, RangeType, SourceType > >(
          over_integrate, diffusion, grid_view, range, source);
}

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion_factor and diffusion_tensor.
 */
template< class DiffusionFactorType, class DiffusionTensorType, class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                           , std::unique_ptr< EllipticLocalizableProduct< DiffusionFactorType, DiffusionTensorType,
                                                                          GridViewType, RangeType, SourceType > >
                           >::type
make_elliptic_localizable_product(const DiffusionFactorType& diffusion_factor,
                                  const DiffusionTensorType& diffusion_tensor,
                                  const GridViewType& grid_view,
                                  const RangeType& range,
                                  const SourceType& source,
                                  const size_t over_integrate = 0)
{
  return DSC::make_unique
      < EllipticLocalizableProduct< DiffusionFactorType, DiffusionTensorType, GridViewType, RangeType, SourceType > >(
          over_integrate, diffusion_factor, diffusion_tensor, grid_view, range, source);
}


// ////////////////////// //
// EllipticMatrixOperator //
// ////////////////////// //

template< class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class RangeSpace,
          class Matrix = typename Stuff::LA::Container< typename RangeSpace::RangeFieldType >::MatrixType,
          class GridView = typename RangeSpace::GridViewType,
          class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType >
class EllipticMatrixOperator
  : public MatrixOperatorDefault< Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume >
{
  typedef MatrixOperatorDefault< Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume > BaseType;
  typedef LocalVolumeIntegralOperator
      < LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > > LocalEllipticOperatorType;
public:
  // see the ctors of EllipticLocalizableProduct
  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value) //
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticMatrixOperator(const DiffusionImp& diffusion, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticMatrixOperator(const size_t over_integrate, const DiffusionImp& diffusion, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticMatrixOperator(const DiffusionFactorImp& diffusion_factor,
                                  const DiffusionTensorImp& diffusion_tensor,
                                  Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticMatrixOperator(const size_t over_integrate,
                                  const DiffusionFactorImp& diffusion_factor,
                                  const DiffusionTensorImp& diffusion_tensor,
                                  Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

private:
  const LocalEllipticOperatorType local_elliptic_operator_;
}; // class EllipticMatrixOperator


// ///////////////////////////// //
// make_elliptic_matrix_operator //
// ///////////////////////////// //

// both diffusion factor and tensor, without matrix

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, space);
\endcode
 */
template< class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      SpaceType, MatrixType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              const SpaceType& space,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType > >(
      over_integrate, diffusion_factor, diffusion_tensor, space);
}

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, space, grid_view);
\endcode
 */
template< class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      SpaceType, MatrixType, GridViewType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              const SpaceType& space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType, GridViewType > >
          (over_integrate, diffusion_factor, diffusion_tensor, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, range_space, source_space, grid_view);
\endcode
 */
template< class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      RangeSpaceType, MatrixType, GridViewType,
                                                                      SourceSpaceType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionFactorType, DiffusionTensorType, RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
          (over_integrate, diffusion_factor, diffusion_tensor, range_space, source_space, grid_view);
}

// both diffusion factor and tensor, with matrix

/**
 * \brief Creates an elliptic matrix operator (source and range space are given by space, grid_view of the space is
 *        used).
 */
template< class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      SpaceType, MatrixType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType > >(
      over_integrate, diffusion_factor, diffusion_tensor, matrix, space);
}

/**
 * \brief Creates an elliptic matrix operator (source and range space are given by space).
 */
template< class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      SpaceType, MatrixType, GridViewType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType, GridViewType > >
          (over_integrate, diffusion_factor, diffusion_tensor, matrix, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix operator.
 */
template< class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionFactorType, DiffusionTensorType,
                                                                      RangeSpaceType, MatrixType, GridViewType,
                                                                      SourceSpaceType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              MatrixType& matrix,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionFactorType, DiffusionTensorType, RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
          (over_integrate, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_view);
}

// single diffusion, without matrix

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically, source and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, space);
\endcode
 */
template< class MatrixType, class DiffusionType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              const SpaceType& space,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType > >
      (over_integrate, diffusion, space);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically, source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, space, grid_view);
\endcode
 */
template< class MatrixType, class DiffusionType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType,
                                                                      GridViewType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              const SpaceType& space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType, GridViewType > >
      (over_integrate, diffusion, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, range_space, source_space, grid_view);
\endcode
 */
template< class MatrixType, class DiffusionType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_matrix< MatrixType >::value
                             && Stuff::is_localizable_function< DiffusionType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, RangeSpaceType, MatrixType,
                                                                      GridViewType, SourceSpaceType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionType, void, RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
          (over_integrate, diffusion, range_space, source_space, grid_view);
}

// single diffusion, with matrix

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, source and range space are given by space,
 *        grid_view of the space is used).
 */
template< class DiffusionType, class MatrixType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType > >(
      over_integrate, diffusion, matrix, space);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, source and range space are given by space).
 */
template< class DiffusionType, class MatrixType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, SpaceType, MatrixType,
                                              GridViewType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionType, void, SpaceType, MatrixType, GridViewType > >
          (over_integrate, diffusion, matrix, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given).
 */
template< class DiffusionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionType >::value
                             && Stuff::LA::is_matrix< MatrixType >::value
                             && is_space< RangeSpaceType >::value
                             && is_space< SourceSpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< EllipticMatrixOperator< DiffusionType, void, RangeSpaceType, MatrixType,
                                                                      GridViewType, SourceSpaceType > >
                           >::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              MatrixType& matrix,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridViewType& grid_view,
                              const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticMatrixOperator
      < DiffusionType, void, RangeSpaceType, MatrixType, GridViewType, SourceSpaceType > >
          (over_integrate, diffusion, matrix, range_space, source_space, grid_view);
}


// //////////////// //
// EllipticOperator //
// //////////////// //

// forward, needed for the traits
template< class DiffusionFactorType, typename DiffusionTensorType, class GridView, class Field = typename DiffusionFactorType::RangeFieldType >
class EllipticOperator;


namespace internal {


template< class DiffusionFactorType, typename DiffusionTensorType, class GridViewType, class Field >
class EllipticOperatorTraits
{
public:
  typedef EllipticOperator< DiffusionFactorType, DiffusionTensorType, GridViewType, Field > derived_type;
  typedef Field FieldType;
};


} // namespace internal


template< class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class GridViewType,
          class Field >
class EllipticOperator
  : public OperatorInterface< internal::EllipticOperatorTraits< DiffusionFactorType, DiffusionTensorType, GridViewType, Field > >
{
  typedef OperatorInterface
      < internal::EllipticOperatorTraits< DiffusionFactorType, DiffusionTensorType, GridViewType, Field > > BaseType;
  typedef LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > LocalEvaluationType;
public:
  using typename BaseType::FieldType;

  template< class ...Args >
  EllipticOperator(const size_t over_integrate, GridViewType grid_view, Args&& ...args)
    : data_functions_(std::forward< Args >(args)...)
    , grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {}

  template< class ...Args >
  EllipticOperator(GridViewType grid_view, Args&& ...args)
    : data_functions_(std::forward< Args >(args)...)
    , grid_view_(grid_view)
    , over_integrate_(0)
  {}

  template< class SourceSpaceType, class VectorType, class RangeSpaceType >
  void apply(const DiscreteFunction< SourceSpaceType, VectorType >& source,
             DiscreteFunction< RangeSpaceType, VectorType >& range) const
  {
    typedef typename Stuff::LA::Container< typename VectorType::ScalarType, VectorType::sparse_matrix_type >::MatrixType
        MatrixType;
    auto op = make_elliptic_matrix_operator< MatrixType >(data_functions_.diffusion_factor(),
                                                          data_functions_.diffusion_tensor(),
                                                          source.space(),
                                                          range.space(),
                                                          grid_view_,
                                                          over_integrate_);
    op->apply(source, range);
  }

  template< class E, class D, size_t d, class R, size_t r, size_t rC >
  FieldType apply2(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& range,
                   const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source) const
  {
    auto product = make_elliptic_localizable_product(data_functions_.diffusion_factor(),
                                                     data_functions_.diffusion_tensor(),
                                                     grid_view_, range, source, over_integrate_);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template< class RangeType, class SourceType >
  void apply_inverse(const RangeType& /*range*/,
                     SourceType& /*source*/,
                     const Stuff::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

  std::vector< std::string > invert_options() const
  {
    DUNE_THROW(NotImplemented, "yet");
    return {"depends_on_the_vector_type_of_the_discrete_function"};
  }

  Stuff::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

private:
  const LocalEvaluationType data_functions_; // We use the local evaluation to store the data functions since it can
  GridViewType grid_view_;                   // handle the case of single diffusion factor, single diffusion tensor and
  const size_t over_integrate_;              // both factor and tensor and creates the required missing data function.
}; // class EllipticOperator


// ////////////////////// //
// make_elliptic_operator //
// ////////////////////// //

template< class GridViewType, class DiffusionType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< DiffusionType >::value
                           , std::unique_ptr< EllipticOperator< DiffusionType, void, GridViewType,
                                                                typename DiffusionType::RangeFieldType > >
                           >::type
make_elliptic_operator(const GridViewType& grid_view,
                       const DiffusionType& diffusion,
                       const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticOperator< DiffusionType, void, GridViewType, typename DiffusionType::RangeFieldType > >
      (over_integrate, grid_view, diffusion);
}

template< class GridViewType, class DiffusionFactorType, class DiffusionTensorType >
    typename std::enable_if<    Stuff::Grid::is_grid_layer< GridViewType >::value
                             && Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                           , std::unique_ptr< EllipticOperator< DiffusionFactorType, DiffusionTensorType, GridViewType,
                                                                typename DiffusionFactorType::RangeFieldType > >
                           >::type
make_elliptic_operator(const GridViewType& grid_view,
                       const DiffusionFactorType& diffusion_factor,
                       const DiffusionTensorType& diffusion_tensor,
                       const size_t over_integrate = 0)
{
  return DSC::make_unique< EllipticOperator< DiffusionFactorType, DiffusionTensorType, GridViewType, typename DiffusionFactorType::RangeFieldType > >
      (over_integrate, grid_view, diffusion_factor, diffusion_tensor);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_HH
