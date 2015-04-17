// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_L2_HH
#define DUNE_GDT_FUNCTIONALS_L2_HH

#include <dune/stuff/common/memory.hh>
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


// forwards
template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType,
          class LocalEvaluationType = LocalEvaluation::Product<FunctionType> >
class L2Volume;

template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp = typename SpaceImp::GridViewType >
class L2Face;


namespace internal {


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp, class LocalEvaluationType >
class L2VolumeTraits
{
  static_assert(Stuff::is_localizable_function< FunctionType >::value,
                "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space< SpaceImp >::value, "SpaceImp has to be derived from SpaceInterface!");
public:
  typedef L2Volume< FunctionType, VectorImp, SpaceImp, GridViewImp, LocalEvaluationType > derived_type;
  typedef VectorImp                                                                       VectorType;
  typedef SpaceImp                                                                        SpaceType;
  typedef GridViewImp                                                                     GridViewType;
  typedef typename VectorType::ScalarType                                                 ScalarType;
}; // class L2VolumeTraits


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2FaceTraits
{
  static_assert(Stuff::is_localizable_function< FunctionType >::value,
                "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space< SpaceImp >::value, "SpaceImp has to be derived from SpaceInterface!");
public:
  typedef L2Face< FunctionType, VectorImp, SpaceImp, GridViewImp > derived_type;
  typedef VectorImp                                                VectorType;
  typedef SpaceImp                                                 SpaceType;
  typedef GridViewImp                                              GridViewType;
  typedef typename VectorType::ScalarType                          ScalarType;
}; // class L2FaceTraits


} // namespace internal


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp, class LocalEvaluationType >
class L2Volume
  : public Functionals::VectorBased< internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp,
                                                               LocalEvaluationType> >
  , public SystemAssembler< SpaceImp, GridViewImp, SpaceImp >
{
  typedef Functionals::VectorBased< internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp,
                                                              LocalEvaluationType> > FunctionalBaseType;
  typedef SystemAssembler< SpaceImp, GridViewImp, SpaceImp >                         AssemblerBaseType;
  typedef LocalFunctional::Codim0Integral<LocalEvaluationType>                       LocalFunctionalType;
  typedef LocalAssembler::Codim0Vector< LocalFunctionalType >                        LocalAssemblerType;

public:
  typedef internal::L2VolumeTraits< FunctionType, VectorImp, SpaceImp, GridViewImp, LocalEvaluationType > Traits;
  typedef typename Traits::VectorType   VectorType;
  typedef typename Traits::SpaceType    SpaceType;
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

  L2Volume(const FunctionType& function, VectorType& vec, const SpaceType& spc,
           const LocalFunctionalType localFunctional)
    : FunctionalBaseType(vec, spc)
    , AssemblerBaseType(spc)
    , function_(function)
    , local_functional_(localFunctional)
    , local_assembler_(local_functional_)
  {
    setup();
  }
  virtual void assemble() override final
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


template< class FunctionType, class VectorImp, class SpaceImp, class GridViewImp >
class L2Face
  : public Functionals::VectorBased< internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
  , public SystemAssembler< SpaceImp, GridViewImp, SpaceImp >
{
  typedef Functionals::VectorBased< internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp > >
                                                                                      FunctionalBaseType;
  typedef SystemAssembler< SpaceImp, GridViewImp, SpaceImp >                          AssemblerBaseType;
  typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > LocalFunctionalType;
  typedef LocalAssembler::Codim1Vector< LocalFunctionalType >                         LocalAssemblerType;
public:
  typedef internal::L2FaceTraits< FunctionType, VectorImp, SpaceImp, GridViewImp >    Traits;

  typedef typename Traits::VectorType   VectorType;
  typedef typename Traits::SpaceType    SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  L2Face(const FunctionType& function,
         VectorType& vec,
         const SpaceType& spc,
         const GridViewType& grd_vw,
         const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* which_intersections
            = new Stuff::Grid::ApplyOn::AllIntersections< GridViewType >())
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
         const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* which_intersections
            = new Stuff::Grid::ApplyOn::AllIntersections< GridViewType >())
    : FunctionalBaseType(vec, spc)
    , AssemblerBaseType(spc)
    , function_(function)
    , local_functional_(function_)
    , local_assembler_(local_functional_)
  {
    setup(which_intersections);
  }

  virtual void assemble() override final
  {
    AssemblerBaseType::assemble();
  }

private:
  void setup(const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* which_intersections)
  {
    this->add(local_assembler_, this->vector(), which_intersections);
  }

  const FunctionType& function_;
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class L2Face


template< class F, class V, class S, class GV >
  std::unique_ptr< L2Volume< F, V, S, GV > >
make_l2_volume(const F& function, V& vector, const S& space, const GV& grid_view)
{
  return Stuff::Common::make_unique< L2Volume< F, V, S, GV > >(function, vector, space, grid_view);
}

template< class F, class V, class S >
  std::unique_ptr< L2Volume< F, V, S > >
make_l2_volume(const F& function, V& vector, const S& space)
{
  return Stuff::Common::make_unique< L2Volume< F, V, S > >(function, vector, space);
}


template< class F, class V, class S >
  std::unique_ptr< L2Face< F, V, S > >
make_l2_face(const F& function,
             V& vector,
             const S& space,
             const Stuff::Grid::ApplyOn::WhichIntersection< typename S::GridViewType >* which_intersections
                = new Stuff::Grid::ApplyOn::AllIntersections< typename S::GridViewType >())
{
  return Stuff::Common::make_unique< L2Face< F, V, S > >(function, vector, space, which_intersections);
}

} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_L2_HH
