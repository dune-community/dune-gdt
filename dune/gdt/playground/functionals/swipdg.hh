// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_FUNCTIONALS_SWIPDG_HH
#define DUNE_GDT_PLAYGROUND_FUNCTIONALS_SWIPDG_HH


//#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/functionals/swipdg.hh>


namespace Dune {
namespace GDT {
namespace Functionals {


template< class DiffusionFactorType
        , class DirichletType
        , class VectorImp
        , class SpaceImp
        , class GridViewImp
        , class DiffusionTensorType >
class DirichletBoundarySWIPDG
  : public Functionals::VectorBased< internal::DirichletBoundarySWIPDGTraits< DiffusionFactorType, DirichletType
                                                                            , VectorImp, SpaceImp, GridViewImp
                                                                            , DiffusionTensorType > >
  , public SystemAssembler< SpaceImp, GridViewImp, SpaceImp >
{
  typedef Functionals::VectorBased< internal::DirichletBoundarySWIPDGTraits< DiffusionFactorType, DirichletType
                                                                           , VectorImp, SpaceImp, GridViewImp
                                                                           , DiffusionTensorType > > FunctionalBaseType;
  typedef SystemAssembler< SpaceImp, GridViewImp, SpaceImp > AssemblerBaseType;

  typedef LocalFunctional::Codim1Integral< LocalEvaluation::SWIPDG::BoundaryRHS< DirichletType, DiffusionFactorType
                                                                               , DiffusionTensorType > >
      LocalFunctionalType;
  typedef LocalAssembler::Codim1Vector< LocalFunctionalType >
      LocalAssemblerType;

  typedef typename VectorImp::ScalarType ScalarType;

public:
  typedef internal::DirichletBoundarySWIPDGTraits
      < DiffusionFactorType, DirichletType, VectorImp, SpaceImp, GridViewImp, DiffusionTensorType > Traits;

  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::SpaceType  SpaceType;
  typedef typename Traits::GridViewType GridViewType;

  typedef Stuff::Grid::BoundaryInfoInterface< typename GridViewType::Intersection > BoundaryInfoType;

  DirichletBoundarySWIPDG(const DiffusionFactorType& diffusion_factor,
                          const DiffusionTensorType& diffusion_tensor,
                          const DirichletType& dirichlet,
                          const BoundaryInfoType& boundary_info,
                          VectorType& vector,
                          const SpaceType& space,
                          const GridViewType& grid_view,
                          const ScalarType beta
                              = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : FunctionalBaseType(vector, space, grid_view)
    , AssemblerBaseType(space, grid_view)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , dirichlet_(dirichlet)
    , boundary_info_(boundary_info)
    , local_functional_(dirichlet_, diffusion_factor_, diffusion_tensor_, beta)
    , local_assembler_(local_functional_)
  {
    setup();
  }

  DirichletBoundarySWIPDG(const DiffusionFactorType& diffusion_factor,
                          const DiffusionTensorType& diffusion_tensor,
                          const DirichletType& dirichlet,
                          const BoundaryInfoType& boundary_info,
                          VectorType& vector,
                          const SpaceType& space,
                          const ScalarType beta
                              = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : FunctionalBaseType(vector, space)
    , AssemblerBaseType(space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , dirichlet_(dirichlet)
    , boundary_info_(boundary_info)
    , local_functional_(dirichlet_, diffusion_factor_, diffusion_tensor_, beta)
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
    this->add(local_assembler_, this->vector(), new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(boundary_info_));
  }

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DirichletType& dirichlet_;
  const BoundaryInfoType& boundary_info_;
  const LocalFunctionalType local_functional_;
  const LocalAssemblerType local_assembler_;
}; // class DirichletBoundarySWIPDG


} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_FUNCTIONALS_SWIPDG_HH
