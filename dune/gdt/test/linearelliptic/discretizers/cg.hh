// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_CG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_CG_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretizations/default.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/projections/dirichlet.hh>
#include <dune/gdt/operators/elliptic.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "../problems/interface.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


/**
 * \brief Discretizes a linear elliptic PDE using a continuous Galerkin Finite Element method.
 */
template <class GridType,
          XT::Grid::Layers layer = XT::Grid::Layers::leaf,
          Backends space_backend = default_cg_backend,
          XT::LA::Backends la = XT::LA::default_sparse_backend,
          int pol = 1,
          class RangeFieldType = double,
          size_t dimRange = 1>
class CGDiscretizer
{
public:
  typedef ProblemInterface<typename GridType::template Codim<0>::Entity,
                           typename GridType::ctype,
                           GridType::dimension,
                           RangeFieldType,
                           dimRange>
      ProblemType;
  typedef CgSpaceProvider<GridType, layer, space_backend, pol, RangeFieldType, dimRange> SpaceProvider;
  typedef typename SpaceProvider::Type SpaceType;
  typedef typename XT::LA::Container<RangeFieldType, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<RangeFieldType, la>::VectorType VectorType;
  typedef StationaryContainerBasedDefaultDiscretization<ProblemType, SpaceType, MatrixType, VectorType, SpaceType>
      DiscretizationType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::cg;
  static const constexpr XT::LA::Backends la_backend = la;
  static const int polOrder = pol;

  static std::string static_id() // int() needed, otherwise we get a
  { // linker error
    return std::string("gdt.linearelliptic.discretization.cg.order_") + Dune::XT::Common::to_string(int(polOrder));
  }

  static DiscretizationType
  discretize(XT::Grid::GridProvider<GridType>& grid_provider, const ProblemType& problem, const int level = 0)
  {
    auto logger = XT::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space = SpaceProvider::create(grid_provider, level);
    logger.debug() << "grid has " << space.grid_layer().indexSet().size(0) << " elements" << std::endl;
    typedef typename SpaceType::GridLayerType GridLayerType;
    typedef typename GridLayerType::Intersection IntersectionType;
    auto boundary_info = XT::Grid::BoundaryInfoFactory<IntersectionType>::create(problem.boundary_info_cfg());
    logger.info() << "Assembling... " << std::endl;
    VectorType rhs_vector(space.mapper().size(), 0.0);
    auto elliptic_operator =
        make_elliptic_matrix_operator<MatrixType>(problem.diffusion_factor(), problem.diffusion_tensor(), space);
    auto l2_force_functional = make_l2_volume_vector_functional(problem.force(), rhs_vector, space);
    auto l2_neumann_functional =
        make_l2_face_vector_functional(problem.neumann(),
                                       rhs_vector,
                                       space,
                                       new XT::Grid::ApplyOn::NeumannIntersections<GridLayerType>(*boundary_info));
    // prepare the dirichlet projection
    auto dirichlet_projection = make_discrete_function<VectorType>(space, "dirichlet values");
    auto dirichlet_projection_operator = make_localizable_dirichlet_projection_operator(
        space.grid_layer(), *boundary_info, problem.dirichlet(), dirichlet_projection);
    DirichletConstraints<IntersectionType> dirichlet_constraints(*boundary_info, space.mapper().size());
    // register everything for assembly in one grid walk
    SystemAssembler<SpaceType> assembler(space);
    assembler.append(*elliptic_operator);
    assembler.append(*l2_force_functional);
    assembler.append(*l2_neumann_functional);
    assembler.append(*dirichlet_projection_operator);
    assembler.append(dirichlet_constraints);
    assembler.assemble();
    // apply the dirichlet shift
    auto& system_matrix = elliptic_operator->matrix();
    rhs_vector -= system_matrix * dirichlet_projection.vector();
    dirichlet_constraints.apply(system_matrix, rhs_vector);
    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, system_matrix, rhs_vector, dirichlet_projection.vector());
  } // ... discretize(...)
}; // class CGDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_CG_HH
