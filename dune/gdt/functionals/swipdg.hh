// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_SWIPDG_HH
#define DUNE_GDT_FUNCTIONALS_SWIPDG_HH

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Functionals {


template <class DiffusionFactorType, class DirichletType, class VectorImp, class SpaceImp,
          class GridViewImp = typename SpaceImp::GridViewType, class DiffusionTensorType = void>
class DirichletBoundarySWIPDG;


namespace internal {


template <class DiffusionFactorType, class DirichletType, class VectorImp, class SpaceImp, class GridViewImp,
          class DiffusionTensorType>
class DirichletBoundarySWIPDGTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<DirichletType>::value,
                "DirichletType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space<SpaceImp>::value, "SpaceImp has to be derived from SpaceInterface!");

public:
  typedef DirichletBoundarySWIPDG<DiffusionFactorType, DirichletType, VectorImp, SpaceImp, GridViewImp,
                                  DiffusionTensorType> derived_type;
  typedef VectorImp VectorType;
  typedef SpaceImp SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
}; // class DirichletBoundarySWIPDGTraits


template <class DiffusionType, class DirichletType, class VectorImp, class SpaceImp, class GridViewImp>
class DirichletBoundarySWIPDGTraits<DiffusionType, DirichletType, VectorImp, SpaceImp, GridViewImp, void>
{
  static_assert(Stuff::is_localizable_function<DiffusionType>::value,
                "DiffusionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<DirichletType>::value,
                "DirichletType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(is_space<SpaceImp>::value, "SpaceImp has to be derived from SpaceInterface!");

public:
  typedef DirichletBoundarySWIPDG<DiffusionType, DirichletType, VectorImp, SpaceImp, GridViewImp, void> derived_type;
  typedef VectorImp VectorType;
  typedef SpaceImp SpaceType;
  typedef GridViewImp GridViewType;
  typedef typename VectorType::ScalarType ScalarType;
}; // class DirichletBoundarySWIPDGTraits< ..., void >


} // namespace internal


template <class DiffusionType, class DirichletType, class VectorImp, class SpaceImp, class GridViewImp>
class DirichletBoundarySWIPDG<DiffusionType, DirichletType, VectorImp, SpaceImp, GridViewImp, void>
    : public Functionals::VectorBased<internal::DirichletBoundarySWIPDGTraits<DiffusionType, DirichletType, VectorImp,
                                                                              SpaceImp, GridViewImp, void>>,
      public SystemAssembler<SpaceImp, GridViewImp, SpaceImp>
{
  typedef Functionals::VectorBased<internal::DirichletBoundarySWIPDGTraits<DiffusionType, DirichletType, VectorImp,
                                                                           SpaceImp, GridViewImp, void>>
      FunctionalBaseType;
  typedef SystemAssembler<SpaceImp, GridViewImp, SpaceImp> AssemblerBaseType;
  typedef LocalFunctional::Codim1Integral<LocalEvaluation::SWIPDG::BoundaryRHS<DirichletType, DiffusionType>>
      LocalFunctionalType;
  typedef LocalAssembler::Codim1Vector<LocalFunctionalType> LocalAssemblerType;
  typedef typename VectorImp::ScalarType ScalarType;

public:
  typedef internal::DirichletBoundarySWIPDGTraits<DiffusionType, DirichletType, VectorImp, SpaceImp, GridViewImp, void>
      Traits;

  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  typedef Stuff::Grid::BoundaryInfoInterface<typename GridViewType::Intersection> BoundaryInfoType;

  DirichletBoundarySWIPDG(
      const DiffusionType& diffusion, const DirichletType& dirichlet, const BoundaryInfoType& boundary_info,
      VectorType& vector, const SpaceType& space, const GridViewType& grid_view,
      const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : FunctionalBaseType(vector, space, grid_view)
    , AssemblerBaseType(space, grid_view)
    , diffusion_(diffusion)
    , dirichlet_(dirichlet)
    , boundary_info_(boundary_info)
    , local_functional_(dirichlet_, diffusion_, beta)
    , local_assembler_(local_functional_)
  {
    setup();
  }

  DirichletBoundarySWIPDG(
      const DiffusionType& diffusion, const DirichletType& dirichlet, const BoundaryInfoType& boundary_info,
      VectorType& vector, const SpaceType& space,
      const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : FunctionalBaseType(vector, space)
    , AssemblerBaseType(space)
    , diffusion_(diffusion)
    , dirichlet_(dirichlet)
    , boundary_info_(boundary_info)
    , local_functional_(dirichlet_, diffusion_, beta)
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
    this->add(local_assembler_,
              this->vector(),
              new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info_));
  }

  const DiffusionType& diffusion_;
  const DirichletType& dirichlet_;
  const BoundaryInfoType& boundary_info_;
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class DirichletBoundarySWIPDG


} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_SWIPDG_HH
