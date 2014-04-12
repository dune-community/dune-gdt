// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Functionals {


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType >
class L2Volume;


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2VolumeTraits
{
  static_assert(std::is_base_of< Stuff::LocalizableFunctionInterface< typename FunctionType::EntityType
                                                                    , typename FunctionType::DomainFieldType
                                                                    , FunctionType::dimDomain
                                                                    , typename FunctionType::RangeFieldType
                                                                    , FunctionType::dimRange
                                                                    , FunctionType::dimRangeCols >
                               , FunctionType >::value,
                "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename SpaceImp::Traits >, SpaceImp >::value,
                "SpaceImp has to be derived from SpaceInterface!");
public:
  typedef L2Volume< FunctionType, VectorImp, SpaceImp, GridViewImp > derived_type;
  typedef VectorImp   VectorType;
  typedef SpaceImp    SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
  typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > LocalFunctionalType;
}; // class L2VolumeTraits


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2Volume
  : public Functionals::AssemblableVolumeBase< L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
{
  typedef Functionals::AssemblableVolumeBase< L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > > BaseType;
public:
  typedef L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > Traits;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType  SpaceType;
  typedef typename Traits::GridViewType GridViewType;

private:
  typedef typename Traits::LocalFunctionalType LocalFunctionalType;

public:
  L2Volume(const FunctionType& function, VectorType& vector, const SpaceType& space, const GridViewType& grid_view)
    : BaseType(vector, space, grid_view)
    , local_functional_(function)
  {}

  L2Volume(const FunctionType& function, VectorType& vector, const SpaceType& space)
    : BaseType(vector, space)
    , local_functional_(function)
  {}

  virtual const LocalFunctionalType& local_functional() const DS_OVERRIDE DS_FINAL
  {
    return local_functional_;
  }

private:
  const LocalFunctionalType local_functional_;
}; // class L2Volume



template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType >
class L2Face;


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2FaceTraits
{
  static_assert(std::is_base_of< Stuff::LocalizableFunctionInterface< typename FunctionType::EntityType
                                                                    , typename FunctionType::DomainFieldType
                                                                    , FunctionType::dimDomain
                                                                    , typename FunctionType::RangeFieldType
                                                                    , FunctionType::dimRange
                                                                    , FunctionType::dimRangeCols >
                               , FunctionType >::value,
                "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename SpaceImp::Traits >, SpaceImp >::value,
                "SpaceImp has to be derived from SpaceInterface!");
public:
  typedef L2Face< FunctionType, VectorImp, SpaceImp, GridViewImp > derived_type;
  typedef VectorImp   VectorType;
  typedef SpaceImp    SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
  typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > LocalFunctionalType;
}; // class L2FaceTraits


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2Face
  : public Functionals::AssemblableFaceBase< L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
{
  typedef Functionals::AssemblableFaceBase< L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > > BaseType;
public:
  typedef L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > Traits;

  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType  SpaceType;
  typedef typename Traits::GridViewType GridViewType;

private:
  typedef typename Traits::LocalFunctionalType LocalFunctionalType;

public:
  L2Face(const FunctionType& function, VectorType& vector, const SpaceType& space, const GridViewType& grid_view)
    : BaseType(vector, space, grid_view)
    , local_functional_(function)
  {}

  L2Face(const FunctionType& function, VectorType& vector, const SpaceType& space)
    : BaseType(vector, space)
    , local_functional_(function)
  {}

  virtual const LocalFunctionalType& local_functional() const DS_OVERRIDE DS_FINAL
  {
    return local_functional_;
  }

private:
  const LocalFunctionalType local_functional_;
}; // class L2Face


} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
