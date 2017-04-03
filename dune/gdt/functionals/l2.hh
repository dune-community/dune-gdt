// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {

// //////////////////////// //
// L2VolumeVectorFunctional //
// //////////////////////// //

/**
 * \todo Unit tests are missing for this class.
 * \note We could add additional ctors which accept XT::Grid::ApplyOn::WhichEntity, analogously to
 *       L2FaceVectorFunctional (see below), but we did not need this until now.
 */
template <class FunctionType,
          class Space,
          class Vector = typename XT::LA::Container<typename Space::RangeFieldType>::VectorType,
          class GridLayer = typename Space::GridLayerType,
          class Field = typename Space::RangeFieldType>
class L2VolumeVectorFunctional : public VectorFunctionalBase<Vector, Space, GridLayer, Field>
{
  typedef L2VolumeVectorFunctional<FunctionType, Space, Vector, GridLayer, Field> ThisType;
  typedef VectorFunctionalBase<Vector, Space, GridLayer, Field> BaseType;

public:
  /// \sa VectorFunctionalBase
  L2VolumeVectorFunctional(const ThisType& other) = delete;
  L2VolumeVectorFunctional(ThisType&& source) = delete;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  template <class... Args>
  explicit L2VolumeVectorFunctional(const FunctionType& function, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(function)
  {
    this->append(local_l2_functional_);
  }

  template <class... Args>
  explicit L2VolumeVectorFunctional(const size_t over_integrate, const FunctionType& function, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->append(local_l2_functional_);
  }

private:
  const LocalVolumeIntegralFunctional<LocalProductIntegrand<FunctionType>> local_l2_functional_;
}; // class L2VolumeVectorFunctional


// //////////////////////////////// //
// make_l2_volume_vector_functional //
// //////////////////////////////// //

template <class VectorType, class FunctionType, class SpaceType>
typename std::enable_if<XT::LA::is_vector<VectorType>::value
                            && XT::Functions::is_localizable_function<FunctionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType>>>::type
make_l2_volume_vector_functional(const FunctionType& function, const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType>>(
      over_integrate, function, space);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_volume_vector_functional(const FunctionType& function,
                                     const SpaceType& space,
                                     const GridLayerType& grid_layer,
                                     const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, function, space, grid_layer);
}

template <class FunctionType, class VectorType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<FunctionType>::value
                            && XT::LA::is_vector<VectorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType>>>::type
make_l2_volume_vector_functional(const FunctionType& function,
                                 VectorType& vector,
                                 const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType>>(
      over_integrate, function, vector, space);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_volume_vector_functional(const FunctionType& function,
                                     VectorType& vector,
                                     const SpaceType& space,
                                     const GridLayerType& grid_layer,
                                     const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2VolumeVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, function, vector, space, grid_layer);
}


// ////////////////////// //
// L2FaceVectorFunctional //
// ////////////////////// //

/**
 * \todo Unit tests are missing for this class.
 */
template <class FunctionType,
          class Space,
          class Vector = typename XT::LA::Container<typename Space::RangeFieldType>::VectorType,
          class GridLayer = typename Space::GridLayerType,
          class Field = typename Space::RangeFieldType>
class L2FaceVectorFunctional : public VectorFunctionalBase<Vector, Space, GridLayer, Field>
{
  typedef VectorFunctionalBase<Vector, Space, GridLayer, Field> BaseType;

public:
  using typename BaseType::GridLayerType;

  template <class... Args>
  explicit L2FaceVectorFunctional(const FunctionType& function, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(function)
  {
    this->append(local_l2_functional_);
  }

  template <class... Args>
  explicit L2FaceVectorFunctional(const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* which_intersections,
                                  const FunctionType& function,
                                  Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(function)
  {
    this->append(local_l2_functional_, which_intersections);
  }

  template <class... Args>
  explicit L2FaceVectorFunctional(const size_t over_integrate, const FunctionType& function, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->append(local_l2_functional_);
  }

  template <class... Args>
  explicit L2FaceVectorFunctional(const size_t over_integrate,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* which_intersections,
                                  const FunctionType& function,
                                  Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->append(local_l2_functional_, which_intersections);
  }

private:
  const LocalFaceIntegralFunctional<LocalProductIntegrand<FunctionType>> local_l2_functional_;
}; // class L2FaceVectorFunctional


// ////////////////////////////// //
// make_l2_face_vector_functional //
// ////////////////////////////// //

template <class VectorType, class FunctionType, class SpaceType>
typename std::enable_if<XT::LA::is_vector<VectorType>::value
                            && XT::Functions::is_localizable_function<FunctionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType>>>::type
make_l2_face_vector_functional(const FunctionType& function, const SpaceType& space, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType>>(
      over_integrate, function, space);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   const SpaceType& space,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      where, function, space);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   const SpaceType& space,
                                   const size_t over_integrate,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, where, function, space);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, function, space, grid_layer);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      where, function, space, grid_layer);
}

template <class VectorType, class FunctionType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<FunctionType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const size_t over_integrate,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, where, function, space, grid_layer);
}

template <class FunctionType, class VectorType, class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<FunctionType>::value
                            && XT::LA::is_vector<VectorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType>>>::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType>>(
      over_integrate, function, vector, space);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   VectorType& vector,
                                   const SpaceType& space,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      where, function, vector, space);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   VectorType& vector,
                                   const SpaceType& space,
                                   const size_t over_integrate,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, where, function, vector, space);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   VectorType& vector,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, function, vector, space, grid_layer);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   VectorType& vector,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      where, function, vector, space, grid_layer);
}

template <class FunctionType, class VectorType, class SpaceType, class GridLayerType>
typename std::
    enable_if<XT::Functions::is_localizable_function<FunctionType>::value && XT::LA::is_vector<VectorType>::value
                  && is_space<SpaceType>::value
                  && XT::Grid::is_layer<GridLayerType>::value,
              std::unique_ptr<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>>::type
    make_l2_face_vector_functional(const FunctionType& function,
                                   VectorType& vector,
                                   const SpaceType& space,
                                   const GridLayerType& grid_layer,
                                   const size_t over_integrate,
                                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where)
{
  return Dune::XT::Common::make_unique<L2FaceVectorFunctional<FunctionType, SpaceType, VectorType, GridLayerType>>(
      over_integrate, where, function, vector, space, grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
