// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretizations/default/stationary-containerbased.hh>
#include <dune/gdt/functionals/elliptic-ipdg.hh>
#include <dune/gdt/functionals/l2.hh>
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
template <class GridType,
          XT::Grid::Layers layer = XT::Grid::Layers::leaf,
          Backends spacebackend = default_dg_backend,
          XT::LA::Backends la = XT::LA::default_sparse_backend,
          int pol = 1,
          class RangeFieldType = double,
          size_t dimRange = 1,
          LocalEllipticIpdgIntegrands::Method method = LocalEllipticIpdgIntegrands::default_method>
class IpdgDiscretizer
{
  template <LocalEllipticIpdgIntegrands::Method in, class Anything = void>
  struct translate
  {
  };

  template <class Anything>
  struct translate<LocalEllipticIpdgIntegrands::Method::sipdg, Anything>
  {
    static const constexpr ChooseDiscretizer value = ChooseDiscretizer::sipdg;
  };

  template <class Anything>
  struct translate<LocalEllipticIpdgIntegrands::Method::swipdg, Anything>
  {
    static const constexpr ChooseDiscretizer value = ChooseDiscretizer::swipdg;
  };

  template <class Anything>
  struct translate<LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor, Anything>
  {
    static const constexpr ChooseDiscretizer value = ChooseDiscretizer::swipdg_affine_factor;
  };

  template <class Anything>
  struct translate<LocalEllipticIpdgIntegrands::Method::swipdg_affine_tensor, Anything>
  {
    static const constexpr ChooseDiscretizer value = ChooseDiscretizer::swipdg_affine_tensor;
  };

public:
  typedef ProblemInterface<typename GridType::template Codim<0>::Entity,
                           typename GridType::ctype,
                           GridType::dimension,
                           RangeFieldType,
                           dimRange>
      ProblemType;
  typedef DgSpaceProvider<GridType, layer, spacebackend, pol, RangeFieldType, dimRange> SpaceProvider;
  typedef typename SpaceProvider::Type SpaceType;
  typedef typename XT::LA::Container<RangeFieldType, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<RangeFieldType, la>::VectorType VectorType;
  typedef StationaryContainerBasedDefaultDiscretization<ProblemType, SpaceType, MatrixType, VectorType, SpaceType>
      DiscretizationType;
  static const constexpr ChooseDiscretizer type = translate<method>::value;
  static const constexpr XT::LA::Backends la_backend = la;
  static const int polOrder = pol;

  static std::string static_id() //                                                        int() needed, otherwise we
  { //                                                                                     get a linker error
    return std::string("gdt.linearelliptic.discretization.swipdg.order_") + Dune::XT::Common::to_string(int(polOrder));
  }

  template <class SubdomainGridType = XT::Grid::none_t>
  static DiscretizationType discretize(XT::Grid::GridProvider<GridType, SubdomainGridType>& grid_provider,
                                       const ProblemType& problem,
                                       const int level = 0)
  {
    auto logger = XT::Common::TimedLogger().get(static_id());
    logger.info() << "Creating space... " << std::endl;
    auto space = SpaceProvider::create(grid_provider, level);
    logger.debug() << "grid has " << space.grid_layer().indexSet().size(0) << " elements" << std::endl;
    typedef typename SpaceType::GridLayerType GridLayerType;
    auto boundary_info = XT::Grid::BoundaryInfoFactory<XT::Grid::extract_intersection_t<GridLayerType>>::create(
        problem.boundary_info_cfg());
    logger.info() << "Assembling... " << std::endl;
    VectorType rhs_vector(space.mapper().size(), 0.0);
    auto ipdg_operator = make_elliptic_ipdg_matrix_operator<MatrixType, method>(
        problem.diffusion_factor(), problem.diffusion_tensor(), *boundary_info, space);
    auto ipdg_boundary_functional = make_elliptic_ipdg_dirichlet_vector_functional<method>(
        problem.dirichlet(), problem.diffusion_factor(), problem.diffusion_tensor(), *boundary_info, rhs_vector, space);
    auto l2_force_functional = make_l2_volume_vector_functional(problem.force(), rhs_vector, space);
    auto l2_neumann_functional =
        make_l2_face_vector_functional(problem.neumann(),
                                       rhs_vector,
                                       space,
                                       new XT::Grid::ApplyOn::NeumannIntersections<GridLayerType>(*boundary_info));
    // register everything for assembly in one grid walk

    ipdg_operator->append(*ipdg_boundary_functional);
    ipdg_operator->append(*l2_force_functional);
    ipdg_operator->append(*l2_neumann_functional);
    ipdg_operator->assemble();
    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, ipdg_operator->matrix(), rhs_vector);
  } // ... discretize(...)

    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, ipdg_operator->matrix(), rhs_vector);
  } // ... discretize(...)
}; // class IpdgDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_IPDG_HH
