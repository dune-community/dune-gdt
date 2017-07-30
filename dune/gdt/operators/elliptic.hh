// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_HH

#include <type_traits>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/spaces/interface.hh>

#include "interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {


// ////////////////////////// //
// EllipticLocalizableProduct //
// ////////////////////////// //

template <class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class GridLayer,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType>
class EllipticLocalizableProduct : public LocalizableProductBase<GridLayer, Range, Source, Field>
{
  typedef LocalizableProductBase<GridLayer, Range, Source, Field> BaseType;

public:
  // Usually, we only have to hold the data functions for the local operator and perfect forward the rest of the
  // arguments to BaseType. Here it is a bit more complicated, since DiffusionTensorType might be void (and
  // DiffusionFactorType is the the only diffusion). To handle this case we require the enable_if hacks below to
  // disable half of the ctors if DiffusionTensorType is void. In addition we require each ctor twice, once with
  // over_integrate and once without (since the perfect forwarding does not allow for default arguments).

  template <typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
            ,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticLocalizableProduct(const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(diffusion, BaseType::parameter())
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
            ,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticLocalizableProduct(const size_t over_integrate, const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(over_integrate, diffusion, BaseType::parameter())
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionFactorImp, // This ctor is only enabled if we are given two diffusion data functions
            typename DiffusionTensorImp, //                                                   (factor and tensor).
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticLocalizableProduct(const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(diffusion_factor, diffusion_tensor, BaseType::parameter())
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionFactorImp, // This ctor is only enabled if we are given two diffusion data functions
            typename DiffusionTensorImp, //                                                   (factor and tensor).
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticLocalizableProduct(const size_t over_integrate,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(over_integrate, diffusion_factor, diffusion_tensor, BaseType::parameter())
  {
    this->append(local_elliptic_operator_);
  }

private:
  const LocalVolumeIntegralOperator<LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType>,
                                    typename Range::LocalfunctionType,
                                    typename Source::LocalfunctionType,
                                    Field>
      local_elliptic_operator_;
}; // class EllipticLocalizableProduct


// ///////////////////////////////// //
// make_elliptic_localizable_product //
// ///////////////////////////////// //

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion.
 */
template <class DiffusionType, class GridLayerType, class RangeType, class SourceType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionType>::value && XT::Grid::is_layer<GridLayerType>::value
                  && XT::Functions::is_localizable_function<RangeType>::value
                  && XT::Functions::is_localizable_function<SourceType>::value,
              std::unique_ptr<EllipticLocalizableProduct<DiffusionType, void, GridLayerType, RangeType, SourceType>>>::
        type
        make_elliptic_localizable_product(const DiffusionType& diffusion,
                                          const GridLayerType& grid_layer,
                                          const RangeType& range,
                                          const SourceType& source,
                                          const size_t over_integrate = 0,
                                          const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::
      make_unique<EllipticLocalizableProduct<DiffusionType, void, GridLayerType, RangeType, SourceType>>(
          over_integrate, diffusion, grid_layer, range, source, param);
}

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion_factor and diffusion_tensor.
 */
template <class DiffusionFactorType, class DiffusionTensorType, class GridLayerType, class RangeType, class SourceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<RangeType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value,
                        std::unique_ptr<EllipticLocalizableProduct<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   GridLayerType,
                                                                   RangeType,
                                                                   SourceType>>>::type
make_elliptic_localizable_product(const DiffusionFactorType& diffusion_factor,
                                  const DiffusionTensorType& diffusion_tensor,
                                  const GridLayerType& grid_layer,
                                  const RangeType& range,
                                  const SourceType& source,
                                  const size_t over_integrate = 0,
                                  const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<EllipticLocalizableProduct<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  GridLayerType,
                                                                  RangeType,
                                                                  SourceType>>(
      over_integrate, diffusion_factor, diffusion_tensor, grid_layer, range, source, param);
}


// ////////////////////// //
// EllipticMatrixOperator //
// ////////////////////// //

template <class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class RangeSpace,
          class Matrix = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridLayer = typename RangeSpace::GridLayerType,
          class SourceSpace = RangeSpace,
          class Field = typename Matrix::ScalarType>
class EllipticMatrixOperator
    : public MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::volume>
{
  typedef EllipticMatrixOperator<DiffusionFactorType,
                                 DiffusionTensorType,
                                 RangeSpace,
                                 Matrix,
                                 GridLayer,
                                 SourceSpace,
                                 Field>
      ThisType;
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::volume> BaseType;

public:
  /// \sa MatrixOperatorBase
  EllipticMatrixOperator(const ThisType& other) = delete;
  EllipticMatrixOperator(ThisType& other) = delete; // <- b.c. of the too perfect forwarding ctor
  EllipticMatrixOperator(ThisType&& source) = delete;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  // see the ctors of EllipticLocalizableProduct
  template <typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
            ,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticMatrixOperator(const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(diffusion)
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
            ,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticMatrixOperator(const size_t over_integrate, const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(over_integrate, diffusion)
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionFactorImp // This ctor is only enabled
            ,
            typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
            ,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticMatrixOperator(const DiffusionFactorImp& diffusion_factor,
                                  const DiffusionTensorImp& diffusion_tensor,
                                  Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(diffusion_factor, diffusion_tensor)
  {
    this->append(local_elliptic_operator_);
  }

  template <typename DiffusionFactorImp // This ctor is only enabled
            ,
            typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
            ,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticMatrixOperator(const size_t over_integrate,
                                  const DiffusionFactorImp& diffusion_factor,
                                  const DiffusionTensorImp& diffusion_tensor,
                                  Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_elliptic_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    this->append(local_elliptic_operator_);
  }

private:
  const LocalVolumeIntegralOperator<LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType>,
                                    typename RangeSpace::BaseFunctionSetType,
                                    typename SourceSpace::BaseFunctionSetType,
                                    Field>
      local_elliptic_operator_;
}; // class EllipticMatrixOperator


// ///////////////////////////// //
// make_elliptic_matrix_operator //
// ///////////////////////////// //

// both diffusion factor and tensor, without matrix

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, space);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                  && is_space<SpaceType>::value,
              std::
                  unique_ptr<EllipticMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType>>>::
        type
        make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                      const DiffusionTensorType& diffusion_tensor,
                                      const SpaceType& space,
                                      const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType>>(
          over_integrate, diffusion_factor, diffusion_tensor, space);
}

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically, source
 *        and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, space, grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionFactorType,
                                                               DiffusionTensorType,
                                                               SpaceType,
                                                               MatrixType,
                                                               GridLayerType>>>::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              const SpaceType& space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionFactorType,
                                                              DiffusionTensorType,
                                                              SpaceType,
                                                              MatrixType,
                                                              GridLayerType>>(
      over_integrate, diffusion_factor, diffusion_tensor, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(factor, tensor, range_space, source_space, grid_layer);
\endcode
 */
template <class MatrixType,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class RangeSpaceType,
          class SourceSpaceType,
          class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionFactorType,
                                                               DiffusionTensorType,
                                                               RangeSpaceType,
                                                               MatrixType,
                                                               GridLayerType,
                                                               SourceSpaceType>>>::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionFactorType,
                                                              DiffusionTensorType,
                                                              RangeSpaceType,
                                                              MatrixType,
                                                              GridLayerType,
                                                              SourceSpaceType>>(
      over_integrate, diffusion_factor, diffusion_tensor, range_space, source_space, grid_layer);
}

// both diffusion factor and tensor, with matrix

/**
 * \brief Creates an elliptic matrix operator (source and range space are given by space, grid_layer of the space is
 *        used).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                  && XT::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value,
              std::
                  unique_ptr<EllipticMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType>>>::
        type
        make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                      const DiffusionTensorType& diffusion_tensor,
                                      MatrixType& matrix,
                                      const SpaceType& space,
                                      const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, MatrixType>>(
          over_integrate, diffusion_factor, diffusion_tensor, matrix, space);
}

/**
 * \brief Creates an elliptic matrix operator (source and range space are given by space).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionFactorType,
                                                               DiffusionTensorType,
                                                               SpaceType,
                                                               MatrixType,
                                                               GridLayerType>>>::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionFactorType,
                                                              DiffusionTensorType,
                                                              SpaceType,
                                                              MatrixType,
                                                              GridLayerType>>(
      over_integrate, diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix operator.
 */
template <class DiffusionFactorType,
          class DiffusionTensorType,
          class MatrixType,
          class RangeSpaceType,
          class SourceSpaceType,
          class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionFactorType,
                                                               DiffusionTensorType,
                                                               RangeSpaceType,
                                                               MatrixType,
                                                               GridLayerType,
                                                               SourceSpaceType>>>::type
make_elliptic_matrix_operator(const DiffusionFactorType& diffusion_factor,
                              const DiffusionTensorType& diffusion_tensor,
                              MatrixType& matrix,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionFactorType,
                                                              DiffusionTensorType,
                                                              RangeSpaceType,
                                                              MatrixType,
                                                              GridLayerType,
                                                              SourceSpaceType>>(
      over_integrate, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_layer);
}

// single diffusion, without matrix

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically, source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, space);
\endcode
 */
template <class MatrixType, class DiffusionType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType>>>::type
make_elliptic_matrix_operator(const DiffusionType& diffusion, const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType>>(
      over_integrate, diffusion, space);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically, source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, space, grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<DiffusionType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType, GridLayerType>>>::type
    make_elliptic_matrix_operator(const DiffusionType& diffusion,
                                  const SpaceType& space,
                                  const GridLayerType& grid_layer,
                                  const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType, GridLayerType>>(
          over_integrate, diffusion, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, MatrixType has to be supllied, a matrix is
 *        created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_matrix_operator< MatrixType >(diffusion, range_space, source_space, grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionType,
                                                               void,
                                                               RangeSpaceType,
                                                               MatrixType,
                                                               GridLayerType,
                                                               SourceSpaceType>>>::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionType,
                                                              void,
                                                              RangeSpaceType,
                                                              MatrixType,
                                                              GridLayerType,
                                                              SourceSpaceType>>(
      over_integrate, diffusion, range_space, source_space, grid_layer);
}

// single diffusion, with matrix

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, source and range space are given by space,
 *        grid_layer of the space is used).
 */
template <class DiffusionType, class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType>>>::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              MatrixType& matrix,
                              const SpaceType& space,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType>>(
      over_integrate, diffusion, matrix, space);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given, source and range space are given by space).
 */
template <class DiffusionType, class MatrixType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionType>::value && XT::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType, GridLayerType>>>::type
    make_elliptic_matrix_operator(const DiffusionType& diffusion,
                                  MatrixType& matrix,
                                  const SpaceType& space,
                                  const GridLayerType& grid_layer,
                                  const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticMatrixOperator<DiffusionType, void, SpaceType, MatrixType, GridLayerType>>(
          over_integrate, diffusion, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix operator (single diffusion given).
 */
template <class DiffusionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticMatrixOperator<DiffusionType,
                                                               void,
                                                               RangeSpaceType,
                                                               MatrixType,
                                                               GridLayerType,
                                                               SourceSpaceType>>>::type
make_elliptic_matrix_operator(const DiffusionType& diffusion,
                              MatrixType& matrix,
                              const RangeSpaceType& range_space,
                              const SourceSpaceType& source_space,
                              const GridLayerType& grid_layer,
                              const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticMatrixOperator<DiffusionType,
                                                              void,
                                                              RangeSpaceType,
                                                              MatrixType,
                                                              GridLayerType,
                                                              SourceSpaceType>>(
      over_integrate, diffusion, matrix, range_space, source_space, grid_layer);
}


// //////////////// //
// EllipticOperator //
// //////////////// //

// forward, needed for the traits
template <class DiffusionFactorType,
          typename DiffusionTensorType,
          class GridLayer,
          class Field = typename DiffusionFactorType::RangeFieldType>
class EllipticOperator;


namespace internal {


template <class DiffusionFactorType, typename DiffusionTensorType, class GridLayerType, class Field>
class EllipticOperatorTraits
{
public:
  typedef EllipticOperator<DiffusionFactorType, DiffusionTensorType, GridLayerType, Field> derived_type;
  typedef NoJacobian JacobianType;
  typedef Field FieldType;
};


} // namespace internal


template <class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class GridLayerType,
          class Field>
class EllipticOperator : public OperatorInterface<internal::EllipticOperatorTraits<DiffusionFactorType,
                                                                                   DiffusionTensorType,
                                                                                   GridLayerType,
                                                                                   Field>>
{
  typedef OperatorInterface<internal::
                                EllipticOperatorTraits<DiffusionFactorType, DiffusionTensorType, GridLayerType, Field>>
      BaseType;
  typedef LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType> LocalIntegrandType;

public:
  using typename BaseType::FieldType;

  template <class... Args>
  EllipticOperator(const size_t over_integrate, GridLayerType grid_layer, Args&&... args)
    : data_functions_(std::forward<Args>(args)...)
    , grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  template <class... Args>
  EllipticOperator(GridLayerType grid_layer, Args&&... args)
    : data_functions_(std::forward<Args>(args)...)
    , grid_layer_(grid_layer)
    , over_integrate_(0)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range,
             const XT::Common::Parameter& param = {}) const
  {
    typedef typename XT::LA::Container<typename VectorType::ScalarType,
                                       VectorType::Traits::sparse_matrix_type>::MatrixType MatrixType;
    auto op = make_elliptic_matrix_operator<MatrixType>(data_functions_.diffusion_factor(),
                                                        data_functions_.diffusion_tensor(),
                                                        source.space(),
                                                        range.space(),
                                                        grid_layer_,
                                                        over_integrate_);
    op->apply(source, range, param);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
                   const XT::Common::Parameter& param = {}) const
  {
    auto product = make_elliptic_localizable_product(data_functions_.diffusion_factor(),
                                                     data_functions_.diffusion_tensor(),
                                                     grid_layer_,
                                                     range,
                                                     source,
                                                     over_integrate_,
                                                     param);
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
  const LocalIntegrandType data_functions_; // We use the local evaluation to store the data functions since it can
  GridLayerType grid_layer_; // handle the case of single diffusion factor, single diffusion tensor and
  const size_t over_integrate_; // both factor and tensor and creates the required missing data function.
}; // class EllipticOperator


// ////////////////////// //
// make_elliptic_operator //
// ////////////////////// //

template <class GridLayerType, class DiffusionType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value,
                        std::unique_ptr<EllipticOperator<DiffusionType,
                                                         void,
                                                         GridLayerType,
                                                         typename DiffusionType::RangeFieldType>>>::type
make_elliptic_operator(const GridLayerType& grid_layer, const DiffusionType& diffusion, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticOperator<DiffusionType, void, GridLayerType, typename DiffusionType::RangeFieldType>>(
          over_integrate, grid_layer, diffusion);
}

template <class GridLayerType, class DiffusionFactorType, class DiffusionTensorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value,
                        std::unique_ptr<EllipticOperator<DiffusionFactorType,
                                                         DiffusionTensorType,
                                                         GridLayerType,
                                                         typename DiffusionFactorType::RangeFieldType>>>::type
make_elliptic_operator(const GridLayerType& grid_layer,
                       const DiffusionFactorType& diffusion_factor,
                       const DiffusionTensorType& diffusion_tensor,
                       const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticOperator<DiffusionFactorType,
                                                        DiffusionTensorType,
                                                        GridLayerType,
                                                        typename DiffusionFactorType::RangeFieldType>>(
      over_integrate, grid_layer, diffusion_factor, diffusion_tensor);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_HH
