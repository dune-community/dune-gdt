// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_HH
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/type_traits.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// //////////////////////////////////////// //
// // EllipticIpdgDirichletVectorFunctional //
// //////////////////////////////////////// //

/**
 * \todo Add tests!
 */
template <class DirichletType,
          class DiffusionFactorType,
          class DiffusionTensorType, // <- may be void
          class Space,
          LocalEllipticIpdgIntegrands::Method method = LocalEllipticIpdgIntegrands::default_method,
          class Vector = typename XT::LA::Container<typename Space::RangeFieldType>::VectorType,
          class GridLayer = typename Space::GridLayerType,
          class Field = typename Space::RangeFieldType>
class EllipticIpdgDirichletVectorFunctional : public VectorFunctionalBase<Vector, Space, GridLayer, Field>
{
  typedef VectorFunctionalBase<Vector, Space, GridLayer, Field> BaseType;
  typedef EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                DiffusionFactorType,
                                                DiffusionTensorType,
                                                Space,
                                                method,
                                                Vector,
                                                GridLayer,
                                                Field>
      ThisType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;

  /// \sa VectorFunctionalBase
  EllipticIpdgDirichletVectorFunctional(const ThisType& other) = delete;
  EllipticIpdgDirichletVectorFunctional(ThisType&& source) = delete;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  /// \name Ctors for single diffusion
  /// \sa The ctors of EllipticLocalizableProduct
  /// \{

  template <typename Diffusion,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<Diffusion, DiffusionFactorType>::value)
                                               && sizeof(Diffusion)>::type,
            class... Args>
  explicit EllipticIpdgDirichletVectorFunctional(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                                 const DirichletType& dirichlet,
                                                 const Diffusion& diffusion,
                                                 Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_functional_(dirichlet, diffusion)
  {
    append_local_functional(boundary_info);
  }

  template <typename Diffusion,
            typename = typename std::enable_if<(std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<Diffusion, DiffusionFactorType>::value)
                                               && sizeof(Diffusion)>::type,
            class... Args>
  explicit EllipticIpdgDirichletVectorFunctional(const size_t over_integrate,
                                                 const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                                 const DirichletType& dirichlet,
                                                 const Diffusion& diffusion,
                                                 Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_functional_(over_integrate, dirichlet, diffusion)
  {
    append_local_functional(boundary_info);
  }

  /// \}
  /// \name Ctors for diffusion factor and tensor
  /// \{

  template <typename DiffusionFactor,
            typename DiffusionTensor,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactor, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactor)>::type,
            class... Args>
  explicit EllipticIpdgDirichletVectorFunctional(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                                 const DirichletType& dirichlet,
                                                 const DiffusionFactor& diffusion_factor,
                                                 const DiffusionTensor& diffusion_tensor,
                                                 Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_functional_(dirichlet, diffusion_factor, diffusion_tensor)
  {
    append_local_functional(boundary_info);
  }

  template <typename DiffusionFactor,
            typename DiffusionTensor,
            typename = typename std::enable_if<(!std::is_same<DiffusionTensorType, void>::value)
                                               && (std::is_same<DiffusionFactor, DiffusionFactorType>::value)
                                               && sizeof(DiffusionFactor)>::type,
            class... Args>
  explicit EllipticIpdgDirichletVectorFunctional(const size_t over_integrate,
                                                 const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info,
                                                 const DirichletType& dirichlet,
                                                 const DiffusionFactor& diffusion_factor,
                                                 const DiffusionTensor& diffusion_tensor,
                                                 Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_functional_(over_integrate, dirichlet, diffusion_factor, diffusion_tensor)
  {
    append_local_functional(boundary_info);
  }
  /// \}

private:
  void append_local_functional(const XT::Grid::BoundaryInfo<IntersectionType>& boundary_info)
  {
    this->append(local_functional_, new XT::Grid::ApplyOn::DirichletIntersections<GridLayerType>(boundary_info));
  }

  const LocalFaceIntegralFunctional<LocalEllipticIpdgIntegrands::
                                        BoundaryRHS<DirichletType, DiffusionFactorType, DiffusionTensorType>,
                                    typename Space::BaseFunctionSetType,
                                    IntersectionType,
                                    Field>
      local_functional_;
}; // class EllipticIpdgDirichletVectorFunctional


// ///////////////////////////////////////////////// //
// // make_elliptic_ipdg_dirichlet_vector_functional //
// ///////////////////////////////////////////////// //

// We have variants for:
// -     vector given / vector not given,
// - single diffusion / both diffusion factor and tensor,
// -  grid layer given / grid layer not given.
// Each of these variants/flavors additionally allows to optionally specify the IPDG method.

// no vector given, both diffusion factor and tensor, no grid layer given, method specified

template <class VectorType,
          LocalEllipticIpdgIntegrands::Method method,
          class DirichletType,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class SpaceType>
typename std::enable_if<XT::LA::is_vector<VectorType>::value
                            && XT::Functions::is_localizable_function<DirichletType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                              DiffusionFactorType,
                                                                              DiffusionTensorType,
                                                                              SpaceType,
                                                                              method,
                                                                              VectorType>>>::type
make_elliptic_ipdg_dirichlet_vector_functional(
    const DirichletType& dirichlet,
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                             DiffusionFactorType,
                                                                             DiffusionTensorType,
                                                                             SpaceType,
                                                                             method,
                                                                             VectorType>>(
      over_integrate, boundary_info, dirichlet, diffusion_factor, diffusion_tensor, space);
}

// no vector given, both diffusion factor and tensor, no grid layer given, default method

template <class VectorType, class DirichletType, class DiffusionFactorType, class DiffusionTensorType, class SpaceType>
typename std::
    enable_if<XT::LA::is_vector<VectorType>::value && XT::Functions::is_localizable_function<DirichletType>::value
                  && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                  && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                  && is_space<SpaceType>::value,
              std::unique_ptr<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                    DiffusionFactorType,
                                                                    DiffusionTensorType,
                                                                    SpaceType,
                                                                    LocalEllipticIpdgIntegrands::default_method,
                                                                    VectorType>>>::type
    make_elliptic_ipdg_dirichlet_vector_functional(
        const DirichletType& dirichlet,
        const DiffusionFactorType& diffusion_factor,
        const DiffusionTensorType& diffusion_tensor,
        const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>&
            boundary_info,
        const SpaceType& space,
        const size_t over_integrate = 0)
{
  return make_elliptic_ipdg_dirichlet_vector_functional<VectorType, LocalEllipticIpdgIntegrands::default_method>(
      dirichlet, diffusion_factor, diffusion_tensor, boundary_info, space, over_integrate);
}

// no vector given, both diffusion factor and tensor, grid layer given, method specified

// ... yet to be implemented!

// no vector given, both diffusion factor and tensor, grid layer given, default method

// ... yet to be implemented!

// no vector given, single diffusion, no grid layer given, method specified

template <class VectorType,
          LocalEllipticIpdgIntegrands::Method method,
          class DirichletType,
          class DiffusionType,
          class SpaceType>
typename std::enable_if<XT::LA::is_vector<VectorType>::value
                            && XT::Functions::is_localizable_function<DirichletType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                              DiffusionType,
                                                                              void,
                                                                              SpaceType,
                                                                              method,
                                                                              VectorType>>>::type
make_elliptic_ipdg_dirichlet_vector_functional(
    const DirichletType& dirichlet,
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                             DiffusionType,
                                                                             void,
                                                                             SpaceType,
                                                                             method,
                                                                             VectorType>>(
      over_integrate, boundary_info, dirichlet, diffusion, space);
}

// no vector given, single diffusion, no grid layer given, default method

// ... yet to be implemented!

// no vector given, single diffusion, grid layer given, method specified

// ... yet to be implemented!

// no vector given, single diffusion, grid layer given, default method

// ... yet to be implemented!

// vector given, both diffusion factor and tensor, no grid layer given, method specified

template <LocalEllipticIpdgIntegrands::Method method,
          class DirichletType,
          class DiffusionFactorType,
          class DiffusionTensorType,
          class VectorType,
          class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DirichletType>::value
                            && XT::Functions::is_localizable_function<DiffusionFactorType>::value
                            && XT::Functions::is_localizable_function<DiffusionTensorType>::value
                            && XT::LA::is_vector<VectorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                              DiffusionFactorType,
                                                                              DiffusionTensorType,
                                                                              SpaceType,
                                                                              method,
                                                                              VectorType>>>::type
make_elliptic_ipdg_dirichlet_vector_functional(
    const DirichletType& dirichlet,
    const DiffusionFactorType& diffusion_factor,
    const DiffusionTensorType& diffusion_tensor,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    VectorType& vector,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                             DiffusionFactorType,
                                                                             DiffusionTensorType,
                                                                             SpaceType,
                                                                             method,
                                                                             VectorType>>(
      over_integrate, boundary_info, dirichlet, diffusion_factor, diffusion_tensor, vector, space);
}

// vector given, both diffusion factor and tensor, no grid layer given, default method

// ... yet to be implemented!

// vector given, both diffusion factor and tensor, grid layer given, method specified

// ... yet to be implemented!

// vector given, both diffusion factor and tensor, grid layer given, default method

// ... yet to be implemented!

// vector given, single diffusion, no grid layer given, method specified

template <LocalEllipticIpdgIntegrands::Method method,
          class DirichletType,
          class DiffusionType,
          class VectorType,
          class SpaceType>
typename std::enable_if<XT::Functions::is_localizable_function<DirichletType>::value
                            && XT::Functions::is_localizable_function<DiffusionType>::value
                            && XT::LA::is_vector<VectorType>::value
                            && is_space<SpaceType>::value,
                        std::unique_ptr<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                              DiffusionType,
                                                                              void,
                                                                              SpaceType,
                                                                              method,
                                                                              VectorType>>>::type
make_elliptic_ipdg_dirichlet_vector_functional(
    const DirichletType& dirichlet,
    const DiffusionType& diffusion,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>>& boundary_info,
    VectorType& vector,
    const SpaceType& space,
    const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<EllipticIpdgDirichletVectorFunctional<DirichletType,
                                                                             DiffusionType,
                                                                             void,
                                                                             SpaceType,
                                                                             method,
                                                                             VectorType>>(
      over_integrate, boundary_info, dirichlet, diffusion, vector, space);
}

// vector given, single diffusion, no grid layer given, default method

// ... yet to be implemented!

// vector given, single diffusion, grid layer given, method specified

// ... yet to be implemented!

// vector given, single diffusion, grid layer given, default method

// ... yet to be implemented!


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_HH
