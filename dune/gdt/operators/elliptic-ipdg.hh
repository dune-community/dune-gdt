// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/type_traits.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// ///////////////////////////// //
// // EllipticIpdgMatrixOperator //
// ///////////////////////////// //

template <class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class RangeSpace,
          LocalEllipticIpdgIntegrands::Method method = LocalEllipticIpdgIntegrands::default_method,
          class Matrix = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridLayer = typename RangeSpace::GridLayerType,
          class SourceSpace = RangeSpace,
          class Field = typename Matrix::ScalarType>
class EllipticIpdgMatrixOperator
    : public MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::face_and_volume>
{
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridLayer, SourceSpace, Field, ChoosePattern::face_and_volume>
      BaseType;
  typedef EllipticIpdgMatrixOperator<DiffusionFactorType,
                                     DiffusionTensorType,
                                     RangeSpace,
                                     method,
                                     Matrix,
                                     GridLayer,
                                     SourceSpace,
                                     Field>
      ThisType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  /// \sa MatrixOperatorBase
  EllipticIpdgMatrixOperator(const ThisType& other) = delete;
  EllipticIpdgMatrixOperator(ThisType& other) = delete; // <- b.c. of the too perfect forwarding ctor
  EllipticIpdgMatrixOperator(ThisType&& source) = delete;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  /// \name Ctors for given single diffusion
  /// \sa The Ctors of EllipticLocalizableProduct.
  /// \{

  virtual ~EllipticIpdgMatrixOperator()
  {
  }

  template <typename DiffusionImp,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionImp& diffusion,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(diffusion)
    , local_coupling_operator_(diffusion)
    , local_boundary_operator_(diffusion)
  {
    append_all(boundary_info);
  }

  template <typename DiffusionImp,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const size_t over_integrate,
                                      const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionImp& diffusion,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(over_integrate, diffusion)
    , local_coupling_operator_(over_integrate, diffusion)
    , local_boundary_operator_(over_integrate, diffusion)
  {
    append_all(boundary_info);
  }

  /// \}
  /// \name Ctors for diffusion factor and tensor
  /// \{

  template <typename DiffusionFactorImp,
            typename DiffusionTensorImp,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(diffusion_factor, diffusion_tensor)
    , local_coupling_operator_(diffusion_factor, diffusion_tensor)
    , local_boundary_operator_(diffusion_factor, diffusion_tensor)
  {
    append_all(boundary_info);
  }

  template <typename DiffusionFactorImp,
            typename DiffusionTensorImp,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const size_t over_integrate,
                                      const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(over_integrate, diffusion_factor, diffusion_tensor)
    , local_coupling_operator_(over_integrate, diffusion_factor, diffusion_tensor)
    , local_boundary_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    append_all(boundary_info);
  }

  /// \}

private:
  void append_all(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info)
  {
    this->append(local_volume_operator_);
    this->append(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridLayerType>());
    this->append(local_boundary_operator_, new XT::Grid::ApplyOn::DirichletIntersections<GridLayerType>(boundary_info));
  }

  typedef typename RangeSpace::BaseFunctionSetType RangeBaseType;
  typedef typename SourceSpace::BaseFunctionSetType SourceBaseType;

  const LocalVolumeIntegralOperator<LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType>,
                                    RangeBaseType,
                                    SourceBaseType,
                                    Field>
      local_volume_operator_;
  const LocalCouplingIntegralOperator<LocalEllipticIpdgIntegrands::
                                          Inner<DiffusionFactorType, DiffusionTensorType, method>,
                                      RangeBaseType,
                                      IntersectionType,
                                      SourceBaseType,
                                      RangeBaseType,
                                      SourceBaseType,
                                      Field>
      local_coupling_operator_;
  const LocalBoundaryIntegralOperator<LocalEllipticIpdgIntegrands::
                                          BoundaryLHS<DiffusionFactorType, DiffusionTensorType, method>,
                                      RangeBaseType,
                                      IntersectionType,
                                      SourceBaseType,
                                      Field>
      local_boundary_operator_;
}; // class EllipticIpdgMatrixOperator


// ////////////////////////////////// //
// make_elliptic_ipdg_matrix_operator //
// ////////////////////////////////// //

// both diffusion factor and tensor, without matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used, source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(factor, tensor, boundary_info, space);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied, source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                              tensor,
                                                                                              boundary_info,
                                                                                              space);
\endcode
 */
template <class MatrixType,
          LocalEllipticIpdgIntegrands::Method method,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, method, MatrixType>>(
          over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used, source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(factor, tensor, boundary_info, space, grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied, source and range space are given by space).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                              tensor,
                                                                                              boundary_info,
                                                                                              space,
                                                                                              grid_layer);
\endcode
 */
template <class MatrixType,
          LocalEllipticIpdgIntegrands::Method method,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class SpaceType,
          class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(factor,
                                                         tensor,
                                                         boundary_info,
                                                         range_space,
                                                         source_space,
                                                         grid_layer);
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
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   RangeSpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, range_space, source_space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                              tensor,
                                                                                              boundary_info,
                                                                                              range_space,
                                                                                              source_space,
                                                                                              grid_layer);
\endcode
 */
template <class MatrixType,
          LocalEllipticIpdgIntegrands::Method method,
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
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   RangeSpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, range_space, source_space, grid_layer);
}

// both diffusion factor and tensor, with matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space, grid_layer of the space is used).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied, source and range space are given by
 *        space, grid_layer of the space is used).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                  tensor,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  space);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class MatrixType,
          class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, method, MatrixType>>(
          over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied, source and range space are given by
 *        space).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                  tensor,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  space,
                                                                                  grid_layer);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class MatrixType,
          class SpaceType,
          class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used).
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
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   RangeSpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(factor,
                                                                                  tensor,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  range_space,
                                                                                  source_space,
                                                                                  grid_layer);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method,
          class DiffusionFactorType,
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
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                   DiffusionTensorType,
                                                                   RangeSpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_layer);
}

// single diffusion, without matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used, source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(diffusion, boundary_info, space);
\endcode
 */
template <class MatrixType, class DiffusionType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType>>(
      over_integrate, boundary_info, diffusion, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied, source and range space are given by space, grid_layer of the space is used).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                              boundary_info,
                                                                                              space);
\endcode
 */
template <class MatrixType, LocalEllipticIpdgIntegrands::Method method, class DiffusionType, class SpaceType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<DiffusionType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType>>>::type
    make_elliptic_ipdg_matrix_operator(
        const DiffusionType& diffusion,
        const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>&
            boundary_info,
        const SpaceType& space,
        const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType>>(
      over_integrate, boundary_info, diffusion, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used, source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(diffusion, boundary_info, space, grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied, source and range space are given by space).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                              boundary_info,
                                                                                              space,
                                                                                              grid_layer);
\endcode
 */
template <class MatrixType,
          LocalEllipticIpdgIntegrands::Method method,
          class DiffusionType,
          class SpaceType,
          class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType, GridLayerType>>(
          over_integrate, boundary_info, diffusion, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType>(diffusion,
                                                         boundary_info,
                                                         range_space,
                                                         source_space,
                                                         grid_layer);
\endcode
 */
template <class MatrixType, class DiffusionType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   RangeSpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion, range_space, source_space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<MatrixType, LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                              boundary_info,
                                                                                              range_space,
                                                                                              source_space,
                                                                                              grid_layer);
\endcode
 */
template <class MatrixType,
          LocalEllipticIpdgIntegrands::Method method,
          class DiffusionType,
          class RangeSpaceType,
          class SourceSpaceType,
          class GridLayerType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   RangeSpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion, range_space, source_space, grid_layer);
}

// single diffusion, with matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space, grid_layer of the space is used).
 */
template <class DiffusionType, class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType>>(
      over_integrate, boundary_info, diffusion, matrix, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied, source and range space are given by
 *        space, grid_layer of the space is used).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  space);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method, class DiffusionType, class MatrixType, class SpaceType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionType>::value && XT::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType>>>::type
    make_elliptic_ipdg_matrix_operator(
        const DiffusionType& diffusion,
        const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>&
            boundary_info,
        MatrixType& matrix,
        const SpaceType& space,
        const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType>>(
      over_integrate, boundary_info, diffusion, matrix, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space).
 */
template <class DiffusionType, class MatrixType, class SpaceType, class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType>>(
      over_integrate, boundary_info, diffusion, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied, source and range space are given by
 *        space).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  space,
                                                                                  grid_layer);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method,
          class DiffusionType,
          class MatrixType,
          class SpaceType,
          class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   SpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const SpaceType& space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticIpdgMatrixOperator<DiffusionType, void, SpaceType, method, MatrixType, GridLayerType>>(
          over_integrate, boundary_info, diffusion, matrix, space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used).
 */
template <class DiffusionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   RangeSpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion, matrix, range_space, source_space, grid_layer);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator<LocalEllipticIpdgIntegrands::swipdg>(diffusion,
                                                                                  boundary_info,
                                                                                  matrix,
                                                                                  range_space,
                                                                                  source_space,
                                                                                  grid_layer);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method,
          class DiffusionType,
          class MatrixType,
          class RangeSpaceType,
          class SourceSpaceType,
          class GridLayerType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_matrix<MatrixType>::value
                            && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionType,
                                                                   void,
                                                                   RangeSpaceType,
                                                                   method,
                                                                   MatrixType,
                                                                   GridLayerType,
                                                                   SourceSpaceType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    MatrixType& matrix,
    const RangeSpaceType& range_space,
    const SourceSpaceType& source_space,
    const GridLayerType& grid_layer,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionType,
                                                                  void,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridLayerType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion, matrix, range_space, source_space, grid_layer);
}

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH
