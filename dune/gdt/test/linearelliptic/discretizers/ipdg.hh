// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/integrals.hh>
#include <dune/gdt/localoperator/integrals.hh>
#include <dune/gdt/spaces/dg.hh>

#include "../problems/interface.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


/**
 * \brief Discretizes a linear elliptic PDE using an interior penalty discontinuous Galerkin Finite Element method.
 */
template< class GridType,
          Stuff::Grid::ChooseLayer layer = Stuff::Grid::ChooseLayer::leaf,
          ChooseSpaceBackend spacebackend = Spaces::default_dg_backend,
          Stuff::LA::ChooseBackend la_backend = Stuff::LA::default_sparse_backend,
          int pol = 1,
          class RangeFieldType = double,
          size_t dimRange = 1,
          LocalEvaluation::EllipticIpdg::Method method = LocalEvaluation::EllipticIpdg::Method::swipdg >
class SwipdgDiscretizer
{
public:
  typedef ProblemInterface< typename GridType::template Codim< 0 >::Entity,
                            typename GridType::ctype,
                            GridType::dimension,
                            RangeFieldType,
                            dimRange > ProblemType;
  typedef Spaces::DGProvider< GridType, layer, spacebackend, pol, RangeFieldType, dimRange > SpaceProvider;
  typedef typename SpaceProvider::Type                                                        SpaceType;
  typedef typename Stuff::LA::Container< RangeFieldType, la_backend >::MatrixType             MatrixType;
  typedef typename Stuff::LA::Container< RangeFieldType, la_backend >::VectorType             VectorType;
  typedef Discretizations::StationaryContainerBasedDefault
      < ProblemType, SpaceType, MatrixType, VectorType, SpaceType >                           DiscretizationType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::swipdg;
  static const int polOrder = pol;

  static std::string static_id() //                                                        int() needed, otherwise we
  { //                                                                                     get a linker error
    return std::string("gdt.linearelliptic.discretization.swipdg.order_") + DSC::to_string(int(polOrder));
  }

  static DiscretizationType discretize(Stuff::Grid::ProviderInterface< GridType >& grid_provider,
                                       const ProblemType& problem,
                                       const int level = 0)
  {
    auto logger = Stuff::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space = SpaceProvider::create(grid_provider, level);
    logger.debug() << "grid has " << space.grid_view().indexSet().size(0) << " elements" << std::endl;
    typedef typename SpaceType::GridViewType    GridViewType;
    typedef typename GridViewType::Intersection IntersectionType;
    auto boundary_info = Stuff::Grid::BoundaryInfoProvider< IntersectionType >::create(problem.boundary_info_cfg());
    logger.info() << "Assembling... " << std::endl;
    VectorType rhs_vector(space.mapper().size(), 0.0);
    MatrixType system_matrix(space.mapper().size(),
                             space.mapper().size(),
                             space.compute_face_and_volume_pattern());
    typedef typename ProblemType::DiffusionFactorType DiffusionFactorType;
    typedef typename ProblemType::DiffusionTensorType DiffusionTensorType;
    typedef typename ProblemType::FunctionType FunctionType;
    // volume terms
    // * lhs
    typedef LocalVolumeIntegralOperator< LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > >
        EllipticOperatorType;
    const EllipticOperatorType                                ellipticOperator(problem.diffusion_factor(),
                                                                               problem.diffusion_tensor());
    const LocalVolumeTwoFormAssembler< EllipticOperatorType > diffusionMatrixAssembler(ellipticOperator);
    // * rhs
    typedef LocalVolumeIntegralFunctional< LocalEvaluation::Product< FunctionType > > ForceFunctionalType;
    const ForceFunctionalType                                   forceFunctional(problem.force());
    const LocalVolumeFunctionalAssembler< ForceFunctionalType > forceVectorAssembler(forceFunctional);
    // inner face terms
    typedef LocalCouplingIntegralOperator<
        LocalEvaluation::EllipticIpdg::Inner< DiffusionFactorType, DiffusionTensorType, method > > CouplingOperatorType;
    const CouplingOperatorType                                  couplingOperator(problem.diffusion_factor(),
                                                                                 problem.diffusion_tensor());
    const LocalCouplingTwoFormAssembler< CouplingOperatorType > couplingMatrixAssembler(couplingOperator);
    // dirichlet boundary face terms
    // * lhs
    typedef LocalBoundaryIntegralOperator< LocalEvaluation::EllipticIpdg::BoundaryLHS
        < DiffusionFactorType, DiffusionTensorType, method > > DirichletOperatorType;
    const DirichletOperatorType                                  dirichletOperator(problem.diffusion_factor(),
                                                                                   problem.diffusion_tensor());
    const LocalBoundaryTwoFormAssembler< DirichletOperatorType > dirichletMatrixAssembler(dirichletOperator);
    // * rhs
    typedef LocalFaceIntegralFunctional< LocalEvaluation::EllipticIpdg::BoundaryRHS
        < FunctionType, DiffusionFactorType, DiffusionTensorType, method > > DirichletFunctionalType;
    const DirichletFunctionalType                                 dirichletFunctional(problem.dirichlet(),
                                                                                      problem.diffusion_factor(),
                                                                                      problem.diffusion_tensor());
    const LocalFaceFunctionalAssembler< DirichletFunctionalType > dirichletVectorAssembler(dirichletFunctional);
    // neumann boundary face terms
    // * rhs
    typedef LocalFaceIntegralFunctional< LocalEvaluation::Product< FunctionType > > NeumannFunctionalType;
    const NeumannFunctionalType                                 neumannFunctional(problem.neumann());
    const LocalFaceFunctionalAssembler< NeumannFunctionalType > neumannVectorAssembler(neumannFunctional);
    // do all the work
    typedef SystemAssembler< SpaceType > SystemAssemblerType;
    SystemAssemblerType systemAssembler(space);
    systemAssembler.add(diffusionMatrixAssembler, system_matrix);
    systemAssembler.add(forceVectorAssembler, rhs_vector);
    systemAssembler.add(couplingMatrixAssembler,
                        system_matrix,
                        new Stuff::Grid::ApplyOn::InnerIntersectionsPrimally< GridViewType >());
    systemAssembler.add(dirichletMatrixAssembler,
                        system_matrix,
                        new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(*boundary_info));
    systemAssembler.add(dirichletVectorAssembler,
                        rhs_vector,
                        new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(*boundary_info));
    systemAssembler.add(neumannVectorAssembler,
                        rhs_vector,
                        new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(*boundary_info));
    systemAssembler.assemble();
    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, system_matrix, rhs_vector);
  } // ... discretize(...)
}; // class SwipdgDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH
