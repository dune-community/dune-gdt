// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2014, 2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <type_traits>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/product.hh>

#include "vector-based.hh"

namespace Dune {
namespace GDT {


/**
 * See also VectorBasedFunctional and FunctionalInterface for a description of the tempate arguments.
 *
 * See also VectorBasedFunctional for the meaning of the different ctors. Some of these will not compile (if GV and
 * AssemblyGridView do not coincide), which is intended.
 *
 * \sa VectorBasedFunctional
 * \sa FunctionalInterface
 */
template <class V, class GV, size_t r = 1, size_t rC = 1, class F = double, class AssemblyGridView = GV>
class L2VolumeVectorFunctional : public VectorBasedFunctional<V, GV, r, rC, F, AssemblyGridView>
{
  using ThisType = L2VolumeVectorFunctional<V, GV, r, rC, F, AssemblyGridView>;
  using BaseType = VectorBasedFunctional<V, GV, r, rC, F, AssemblyGridView>;

public:
  using typename BaseType::ApplyOnAllElements;
  using typename BaseType::DofFieldType;
  using typename BaseType::ElementFilterType;
  using typename BaseType::ElementType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceVectorType;

  using InducingFunctionType = XT::Functions::GridFunctionInterface<ElementType, r, rC, F>;

private:
  using LocalFunctionalType = LocalElementIntegralFunctional<ElementType, r, rC, F, DofFieldType>;
  using LocalIntegrandType = LocalElementProductIntegrand<ElementType, r, rC, F, DofFieldType>;

public:
  L2VolumeVectorFunctional(AssemblyGridView assembly_grid_view,
                           const SourceSpaceType& source_spc,
                           SourceVectorType& vec,
                           const InducingFunctionType& inducing_function,
                           const size_t over_integrate = 0,
                           const XT::Common::Parameter& param = {},
                           const ElementFilterType& filter = ApplyOnAllElements())
    : BaseType(assembly_grid_view, source_spc, vec)
  {
    this->append(LocalFunctionalType(local_binary_to_unary_element_integrand(inducing_function, LocalIntegrandType(1)),
                                     over_integrate),
                 param,
                 filter);
  }

  L2VolumeVectorFunctional(AssemblyGridView assembly_grid_view,
                           const SourceSpaceType& source_spc,
                           const InducingFunctionType& inducing_function,
                           const size_t over_integrate = 0,
                           const XT::Common::Parameter& param = {},
                           const ElementFilterType& filter = ApplyOnAllElements())
    : BaseType(assembly_grid_view, source_spc)
  {
    this->append(LocalFunctionalType(local_binary_to_unary_element_integrand(inducing_function, LocalIntegrandType(1)),
                                     over_integrate),
                 param,
                 filter);
  }
}; // class L2VolumeVectorFunctional


/// \name Variants of make_l2_volume_vector_functional which accept an existing vector into which to assemble.
/// \{

template <class AGV, class GV, size_t r, size_t rC, class V, class F>
L2VolumeVectorFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F, GridView<AGV>>
make_l2_volume_vector_functional(
    GridView<AGV> assembly_grid_view,
    const SpaceInterface<GV, r, rC, F>& space,
    XT::LA::VectorInterface<V>& vector,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, F>& inducing_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<GridView<AGV>>& filter = XT::Grid::ApplyOn::AllElements<GridView<AGV>>())
{
  return L2VolumeVectorFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F, GridView<AGV>>(
      assembly_grid_view, space, vector.as_imp(), inducing_function, over_integrate, param, filter);
}

template <class GV, size_t r, size_t rC, class F, class V>
L2VolumeVectorFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>
make_l2_volume_vector_functional(
    const SpaceInterface<GV, r, rC, F>& space,
    XT::LA::VectorInterface<V>& vector,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, F>& inducing_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<GV>& filter = XT::Grid::ApplyOn::AllElements<GV>())
{
  return L2VolumeVectorFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>(
      space.grid_view(), space, vector.as_imp(), inducing_function, over_integrate, param, filter);
}

/// \}
/// \name Variants of make_l2_volume_vector_functional which create an appropriate vector into which to assemble.
/// \{

template <class VectorType, class GV, size_t r, size_t rC, class F, class AGV>
typename std::enable_if<XT::LA::is_vector<VectorType>::value,
                        L2VolumeVectorFunctional<VectorType, GV, r, rC, F, AGV>>::type
make_l2_volume_vector_functional(
    GridView<AGV> assembly_grid_view,
    const SpaceInterface<GV, r, rC, F>& space,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, F>& inducing_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<AGV>& filter = XT::Grid::ApplyOn::AllElements<AGV>())
{
  return L2VolumeVectorFunctional<VectorType, GV, r, rC, F, AGV>(
      assembly_grid_view, space, inducing_function, over_integrate, param, filter);
}

template <class VectorType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_vector<VectorType>::value, L2VolumeVectorFunctional<VectorType, GV, r, rC, F>>::type
make_l2_volume_vector_functional(
    const SpaceInterface<GV, r, rC, F>& space,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GV>, r, rC, F>& inducing_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<GV>& filter = XT::Grid::ApplyOn::AllElements<GV>())
{
  return L2VolumeVectorFunctional<VectorType, GV, r, rC, F>(
      space.grid_view(), space, inducing_function, over_integrate, param, filter);
}

/// \}


#if 0
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
class L2FaceVectorFunctional : public VectorBasedFunctional<Vector, Space, GridLayer, Field>
{
  typedef VectorBasedFunctional<Vector, Space, GridLayer, Field> BaseType;

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
  const LocalFaceIntegralFunctional<LocalProductIntegrand<FunctionType>,
                                    typename Space::BaseFunctionSetType,
                                    XT::Grid::extract_intersection_t<GridLayer>,
                                    Field>
      local_l2_functional_;
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
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
