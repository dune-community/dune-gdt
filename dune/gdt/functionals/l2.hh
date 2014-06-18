// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Functionals {


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType,
          class LocalEvaluationType = LocalEvaluation::Product<FunctionType> >
class L2Volume;


namespace internal {


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp, class LocalEvaluationType >
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
  typedef L2Volume< FunctionType, VectorImp, SpaceImp, GridViewImp, LocalEvaluationType > derived_type;
  typedef VectorImp   VectorType;
  typedef SpaceImp    SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
}; // class L2VolumeTraits


} // namespace internal


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp, class LocalEvaluationType >
class L2Volume
  : public Functionals::VectorBased< internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp,
                                                               LocalEvaluationType> >
  , public SystemAssembler< SpaceImp, GridViewImp, SpaceImp >
{
  typedef Functionals::VectorBased< internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp,
                                                              LocalEvaluationType> >
      FunctionalBaseType;
  typedef SystemAssembler< SpaceImp, GridViewImp, SpaceImp > AssemblerBaseType;

  typedef LocalFunctional::Codim0Integral<LocalEvaluationType> LocalFunctionalType;
  typedef LocalAssembler::Codim0Vector< LocalFunctionalType > LocalAssemblerType;

public:
  typedef internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp, LocalEvaluationType > Traits;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType  SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  L2Volume(const FunctionType& function, VectorType& vec, const SpaceType& spc, const GridViewType& grd_vw)
    : FunctionalBaseType(vec, spc, grd_vw)
    , AssemblerBaseType(spc, grd_vw)
    , function_(function)
    , local_functional_(function_)
    , local_assembler_(local_functional_)
  {
    setup();
  }

  L2Volume(const FunctionType& function, VectorType& vec, const SpaceType& spc)
    : FunctionalBaseType(vec, spc)
    , AssemblerBaseType(spc)
    , function_(function)
    , local_functional_(function_)
    , local_assembler_(local_functional_)
  {
    setup();
  }

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  void setup()
  {
    this->add(local_assembler_, this->vector());
  }

  const FunctionType& function_;
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class L2Volume


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType >
class L2Face;


namespace internal {


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
}; // class L2FaceTraits


} // namespace internal


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2Face
  : public Functionals::VectorBased< internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
  , public SystemAssembler< SpaceImp, GridViewImp, SpaceImp >
{
  typedef Functionals::VectorBased< internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
      FunctionalBaseType;
  typedef SystemAssembler< SpaceImp, GridViewImp, SpaceImp > AssemblerBaseType;

  typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > LocalFunctionalType;
  typedef LocalAssembler::Codim1Vector< LocalFunctionalType > LocalAssemblerType;
public:
  typedef internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > Traits;

  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType  SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  L2Face(const FunctionType& function,
         VectorType& vec,
         const SpaceType& spc,
         const GridViewType& grd_vw,
         const GDT::ApplyOn::WhichIntersection< GridViewType >* which_intersections
            = new GDT::ApplyOn::AllIntersections< GridViewType >())
    : FunctionalBaseType(vec, spc, grd_vw)
    , AssemblerBaseType(spc, grd_vw)
    , function_(function)
    , local_functional_(function_)
    , local_assembler_(local_functional_)
  {
    setup(which_intersections);
  }

  L2Face(const FunctionType& function,
         VectorType& vec,
         const SpaceType& spc,
         const GDT::ApplyOn::WhichIntersection< GridViewType >* which_intersections
            = new GDT::ApplyOn::AllIntersections< GridViewType >())
    : FunctionalBaseType(vec, spc)
    , AssemblerBaseType(spc)
    , function_(function)
    , local_functional_(function_)
    , local_assembler_(local_functional_)
  {
    setup(which_intersections);
  }

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  void setup(const GDT::ApplyOn::WhichIntersection< GridViewType >* which_intersections)
  {
    this->add(local_assembler_, this->vector(), which_intersections);
  }

  const FunctionType& function_;
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class L2Face


} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
