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
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/functionals/elliptic-ipdg.hh>
#include <dune/gdt/operators/elliptic-ipdg.hh>
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
    auto ipdg_operator = make_elliptic_ipdg_matrix_operator< MatrixType, method >(problem.diffusion_factor(),
                                                                                  problem.diffusion_tensor(),
                                                                                  *boundary_info,
                                                                                  space);
    auto ipdg_boundary_functional = make_elliptic_ipdg_dirichlet_vector_functional< method >(problem.dirichlet(),
                                                                                             problem.diffusion_factor(),
                                                                                             problem.diffusion_tensor(),
                                                                                             *boundary_info,
                                                                                             rhs_vector,
                                                                                             space);
    auto l2_force_functional = make_l2_volume_vector_functional(problem.force(), rhs_vector, space);
    auto l2_neumann_functional
        = make_l2_face_vector_functional(problem.neumann(),
                                         rhs_vector,
                                         space,
                                         new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(*boundary_info));
    // register everything for assembly in one grid walk
    SystemAssembler< SpaceType > assembler(space);
    assembler.add(*ipdg_operator);
    assembler.add(*ipdg_boundary_functional);
    assembler.add(*l2_force_functional);
    assembler.add(*l2_neumann_functional);
    assembler.assemble();
    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, ipdg_operator->matrix(), rhs_vector);
  } // ... discretize(...)
}; // class SwipdgDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH
