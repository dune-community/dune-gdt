// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// ///////////////////////////// //
// // EllipticIpdgMatrixOperator //
// ///////////////////////////// //

template <class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class RangeSpace, LocalEllipticIpdgIntegrands::Method method = LocalEllipticIpdgIntegrands::default_method,
          class Matrix   = typename XT::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridView = typename RangeSpace::GridViewType, class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class EllipticIpdgMatrixOperator
    : public MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::face_and_volume>
{
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::face_and_volume> BaseType;
  typedef LocalVolumeIntegralOperator<LocalEllipticIntegrand<DiffusionFactorType, DiffusionTensorType>>
      LocalVolumeOperatorType;
  typedef LocalCouplingIntegralOperator<LocalEllipticIpdgIntegrands::Inner<DiffusionFactorType, DiffusionTensorType,
                                                                           method>>
      LocalCouplingOperatorType;
  typedef LocalBoundaryIntegralOperator<LocalEllipticIpdgIntegrands::BoundaryLHS<DiffusionFactorType,
                                                                                 DiffusionTensorType, method>>
      LocalBoundaryOperatorType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::IntersectionType;

  /// \name Ctors for given single diffusion
  /// \sa The Ctors of EllipticLocalizableProduct.
  /// \{

  template <typename DiffusionImp,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(diffusion)
    , local_coupling_operator_(diffusion)
    , local_boundary_operator_(diffusion)
  {
    this->add(local_volume_operator_);
    this->add(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(local_boundary_operator_, new XT::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info));
  }

  template <typename DiffusionImp,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const size_t over_integrate,
                                      const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionImp& diffusion, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(over_integrate, diffusion)
    , local_coupling_operator_(over_integrate, diffusion)
    , local_boundary_operator_(over_integrate, diffusion)
  {
    this->add(local_volume_operator_);
    this->add(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(local_boundary_operator_, new XT::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info));
  }

  /// \}
  /// \name Ctors for diffusion factor and tensor
  /// \{

  template <typename DiffusionFactorImp, typename DiffusionTensorImp,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(diffusion_factor, diffusion_tensor)
    , local_coupling_operator_(diffusion_factor, diffusion_tensor)
    , local_boundary_operator_(diffusion_factor, diffusion_tensor)
  {
    this->add(local_volume_operator_);
    this->add(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(local_boundary_operator_, new XT::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info));
  }

  template <typename DiffusionFactorImp, typename DiffusionTensorImp,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactorImp, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactorImp)>::type,
            class... Args>
  explicit EllipticIpdgMatrixOperator(const size_t over_integrate,
                                      const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_volume_operator_(over_integrate, diffusion_factor, diffusion_tensor)
    , local_coupling_operator_(over_integrate, diffusion_factor, diffusion_tensor)
    , local_boundary_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    this->add(local_volume_operator_);
    this->add(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(local_boundary_operator_, new XT::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info));
  }

  /// \}

private:
  const LocalVolumeOperatorType local_volume_operator_;
  const LocalCouplingOperatorType local_coupling_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class EllipticIpdgMatrixOperator


// ////////////////////////////////// //
// make_elliptic_ipdg_matrix_operator //
// ////////////////////////////////// //

// both diffusion factor and tensor, without matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used, source and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< MatrixType >(factor, tensor, boundary_info, space);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<typename SpaceType::GridViewType::Intersection>& boundary_info, const SpaceType& space,
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
 *        IPDG method has to be supplied, source and range space are given by space, grid_view of the space is used).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< MatrixType, LocalEllipticIpdgIntegrands::swipdg >(factor, tensor,
boundary_info, space);
\endcode
 */
template <class MatrixType, LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType,
          class DiffusionTensorType, class SpaceType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   method, MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<typename SpaceType::GridViewType::Intersection>& boundary_info, const SpaceType& space,
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
auto op = make_elliptic_ipdg_matrix_operator< MatrixType >(factor, tensor, boundary_info, space, grid_view);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType, class GridViewType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType, GridViewType>>>::type
make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                   const DiffusionTensorType& diffusion_tensor,
                                   const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                   const SpaceType& space, const GridViewType& grid_view,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridViewType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied, source and range space are given by space).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< MatrixType, LocalEllipticIpdgIntegrands::swipdg >(factor, tensor,
boundary_info, space, grid_view);
\endcode
 */
template <class MatrixType, LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType,
          class DiffusionTensorType, class SpaceType, class GridViewType>
typename std::enable_if<XT::LA::is_matrix<MatrixType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   method, MatrixType, GridViewType>>>::type
make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                   const DiffusionTensorType& diffusion_tensor,
                                   const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                   const SpaceType& space, const GridViewType& grid_view,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridViewType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        default IPDG method is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< MatrixType >(factor, tensor, boundary_info, range_space, source_space,
grid_view);
\endcode
 */
template <class MatrixType, class DiffusionFactorType, class DiffusionTensorType, class RangeSpaceType,
          class SourceSpaceType, class GridViewType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, RangeSpaceType,
                                                         LocalEllipticIpdgIntegrands::default_method, MatrixType,
                                                         GridViewType, SourceSpaceType>>>::type
    make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                       const DiffusionTensorType& diffusion_tensor,
                                       const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                       const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                       const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridViewType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, range_space, source_space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (MatrixType has to be supllied, a matrix is created automatically,
 *        IPDG method has to be supplied).
 * \note  MatrixType and IPDG method have to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< MatrixType, LocalEllipticIpdgIntegrands::swipdg >(factor, tensor,
boundary_info, range_space, source_space, grid_view);
\endcode
 */
template <class MatrixType, LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType,
          class DiffusionTensorType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::
    enable_if<XT::LA::is_matrix<MatrixType>::value && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value && is_space<RangeSpaceType>::value
                  && is_space<SourceSpaceType>::value && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, RangeSpaceType,
                                                         method, MatrixType, GridViewType, SourceSpaceType>>>::type
    make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                       const DiffusionTensorType& diffusion_tensor,
                                       const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                       const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                       const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridViewType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, range_space, source_space, grid_view);
}

// both diffusion factor and tensor, with matrix

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space, grid_view of the space is used).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<typename SpaceType::GridViewType::Intersection>& boundary_info, MatrixType& matrix,
    const SpaceType& space, const size_t over_integrate = 0)
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
 *        space, grid_view of the space is used).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< LocalEllipticIpdgIntegrands::swipdg >(factor, tensor, boundary_info,
matrix, space);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType, class DiffusionTensorType,
          class MatrixType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   method, MatrixType>>>::type
make_elliptic_ipdg_matrix_operator(
    const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<typename SpaceType::GridViewType::Intersection>& boundary_info, MatrixType& matrix,
    const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType, method, MatrixType>>(
          over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used, source and range space are given by
 *        space).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class SpaceType, class GridViewType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   LocalEllipticIpdgIntegrands::default_method,
                                                                   MatrixType, GridViewType>>>::type
make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                   const DiffusionTensorType& diffusion_tensor,
                                   const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                   MatrixType& matrix, const SpaceType& space, const GridViewType& grid_view,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridViewType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied, source and range space are given by
 *        space).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< LocalEllipticIpdgIntegrands::swipdg >(factor, tensor, boundary_info,
matrix, space, grid_view);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType, class DiffusionTensorType,
          class MatrixType, class SpaceType, class GridViewType>
typename std::enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value
                            && XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, SpaceType,
                                                                   method, MatrixType, GridViewType>>>::type
make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                   const DiffusionTensorType& diffusion_tensor,
                                   const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                   MatrixType& matrix, const SpaceType& space, const GridViewType& grid_view,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  SpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridViewType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (default IPDG method is used).
 */
template <class DiffusionFactorType, class DiffusionTensorType, class MatrixType, class RangeSpaceType,
          class SourceSpaceType, class GridViewType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value && XT::LA::is_matrix<MatrixType>::value
                  && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
                  && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, RangeSpaceType,
                                                         LocalEllipticIpdgIntegrands::default_method, MatrixType,
                                                         GridViewType, SourceSpaceType>>>::type
    make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                       const DiffusionTensorType& diffusion_tensor,
                                       const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                       MatrixType& matrix, const RangeSpaceType& range_space,
                                       const SourceSpaceType& source_space, const GridViewType& grid_view,
                                       const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  LocalEllipticIpdgIntegrands::default_method,
                                                                  MatrixType,
                                                                  GridViewType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_view);
}

/**
 * \brief Creates an elliptic matrix IPDG operator (IPDG method has to be supplied).
 * \note  IPDG method has to be supplied, i.e., use like
\code
auto op = make_elliptic_ipdg_matrix_operator< LocalEllipticIpdgIntegrands::swipdg >(factor, tensor, boundary_info,
matrix, range_space, source_space, grid_view);
\endcode
 */
template <LocalEllipticIpdgIntegrands::Method method, class DiffusionFactorType, class DiffusionTensorType,
          class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::
    enable_if<XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value && XT::LA::is_matrix<MatrixType>::value
                  && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
                  && XT::Grid::is_layer<GridViewType>::value,
              std::unique_ptr<EllipticIpdgMatrixOperator<DiffusionFactorType, DiffusionTensorType, RangeSpaceType,
                                                         method, MatrixType, GridViewType, SourceSpaceType>>>::type
    make_elliptic_ipdg_matrix_operator(const DiffusionFactorType& diffusion_factor,
                                       const DiffusionTensorType& diffusion_tensor,
                                       const XT::Grid::BoundaryInfo<typename GridViewType::Intersection>& boundary_info,
                                       MatrixType& matrix, const RangeSpaceType& range_space,
                                       const SourceSpaceType& source_space, const GridViewType& grid_view,
                                       const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgMatrixOperator<DiffusionFactorType,
                                                                  DiffusionTensorType,
                                                                  RangeSpaceType,
                                                                  method,
                                                                  MatrixType,
                                                                  GridViewType,
                                                                  SourceSpaceType>>(
      over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, range_space, source_space, grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_OPERATORS_ELLIPTIC_IPDG_HH
