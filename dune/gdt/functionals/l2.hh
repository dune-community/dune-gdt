// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/localfunctional/integrals.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/assembler/system.hh>

#include "default.hh"

namespace Dune {
namespace GDT {

// //////////////////////// //
// L2VolumeVectorFunctional //
// //////////////////////// //

/**
 * \todo Unit tests are missing for this class.
 * \note We could add additional ctors which accept Stuff::Grid::ApplyOn::WhichEntity, analogously to
 *       L2FaceVectorFunctional (see below), but we did not need this until now.
 */
template< class FunctionType,
          class Space,
          class Vector = typename Stuff::LA::Container< typename Space::RangeFieldType >::VectorType,
          class GridView = typename Space::GridViewType,
          class Field = typename Space::RangeFieldType>
class L2VolumeVectorFunctional
  : public VectorFunctionalDefault< Vector, Space, GridView, Field >
{
  typedef VectorFunctionalDefault< Vector, Space, GridView, Field > BaseType;
public:
  template< class ...Args >
  explicit L2VolumeVectorFunctional(const FunctionType& function, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(function)
  {
    this->add(local_l2_functional_);
  }

  template< class ...Args >
  explicit L2VolumeVectorFunctional(const size_t over_integrate, const FunctionType& function, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->add(local_l2_functional_);
  }

private:
  const LocalVolumeIntegralFunctional< LocalEvaluation::Product< FunctionType > > local_l2_functional_;
}; // class L2VolumeVectorFunctional


// //////////////////////////////// //
// make_l2_volume_vector_functional //
// //////////////////////////////// //

template< class VectorType, class FunctionType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_volume_vector_functional(const FunctionType& function, const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, function, space);
}

template< class VectorType, class FunctionType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_volume_vector_functional(const FunctionType& function,
                                 const SpaceType& space,
                                 const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, function, space, grid_view);
}

template< class FunctionType, class VectorType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_volume_vector_functional(const FunctionType& function,
                                 VectorType& vector,
                                 const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, function, vector, space);
}

template< class FunctionType, class VectorType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_volume_vector_functional(const FunctionType& function,
                                 VectorType& vector,
                                 const SpaceType& space,
                                 const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique< L2VolumeVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, function, vector, space, grid_view);
}


// ////////////////////// //
// L2FaceVectorFunctional //
// ////////////////////// //

/**
 * \todo Unit tests are missing for this class.
 */
template< class FunctionType,
          class Space,
          class Vector = typename Stuff::LA::Container< typename Space::RangeFieldType >::VectorType,
          class GridView = typename Space::GridViewType,
          class Field = typename Space::RangeFieldType>
class L2FaceVectorFunctional
  : public VectorFunctionalDefault< Vector, Space, GridView, Field >
{
  typedef VectorFunctionalDefault< Vector, Space, GridView, Field > BaseType;
public:
  using typename BaseType::GridViewType;

  template< class ...Args >
  explicit L2FaceVectorFunctional(const FunctionType& function, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(function)
  {
    this->add(local_l2_functional_);
  }

  template< class ...Args >
  explicit L2FaceVectorFunctional(const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* which_intersections,
                                  const FunctionType& function,
                                  Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(function)
  {
    this->add(local_l2_functional_, which_intersections);
  }

  template< class ...Args >
  explicit L2FaceVectorFunctional(const size_t over_integrate, const FunctionType& function, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->add(local_l2_functional_);
  }

  template< class ...Args >
  explicit L2FaceVectorFunctional(const size_t over_integrate,
                                  const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* which_intersections,
                                  const FunctionType& function,
                                  Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_l2_functional_(over_integrate, function)
  {
    this->add(local_l2_functional_, which_intersections);
  }

private:
  const LocalFaceIntegralFunctional< LocalEvaluation::Product< FunctionType > > local_l2_functional_;
}; // class L2FaceVectorFunctional


// ////////////////////////////// //
// make_l2_face_vector_functional //
// ////////////////////////////// //

template< class VectorType, class FunctionType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function, const SpaceType& space, const size_t over_integrate = 0)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, function, space);
}

template< class VectorType, class FunctionType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               const SpaceType& space,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      where, function, space);
}

template< class VectorType, class FunctionType, class SpaceType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               const SpaceType& space,
                               const size_t over_integrate,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, where, function, space);
}

template< class VectorType, class FunctionType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const size_t over_integrate = 0)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, function, space, grid_view);
}

template< class VectorType, class FunctionType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      where, function, space, grid_view);
}

template< class VectorType, class FunctionType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::LA::is_vector< VectorType >::value
                             && Stuff::is_localizable_function< FunctionType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const size_t over_integrate,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, where, function, space, grid_view);
}

template< class FunctionType, class VectorType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const size_t over_integrate = 0)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, function, vector, space);
}

template< class FunctionType, class VectorType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      where, function, vector, space);
}

template< class FunctionType, class VectorType, class SpaceType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const size_t over_integrate,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType > >(
      over_integrate, where, function, vector, space);
}

template< class FunctionType, class VectorType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const size_t over_integrate = 0)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, function, vector, space, grid_view);
}

template< class FunctionType, class VectorType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      where, function, vector, space, grid_view);
}

template< class FunctionType, class VectorType, class SpaceType, class GridViewType >
    typename std::enable_if<    Stuff::is_localizable_function< FunctionType >::value
                             && Stuff::LA::is_vector< VectorType >::value
                             && is_space< SpaceType >::value
                             && Stuff::Grid::is_grid_layer< GridViewType >::value
                           , std::unique_ptr< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > > >::type
make_l2_face_vector_functional(const FunctionType& function,
                               VectorType& vector,
                               const SpaceType& space,
                               const GridViewType& grid_view,
                               const size_t over_integrate,
                               const Stuff::Grid::ApplyOn::WhichIntersection< typename SpaceType::GridViewType >* where)
{
  return DSC::make_unique< L2FaceVectorFunctional< FunctionType, SpaceType, VectorType, GridViewType > >(
      over_integrate, where, function, vector, space, grid_view);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
