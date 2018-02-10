// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BLOCK_IPDG_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BLOCK_IPDG_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/view/subdomain/view.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/operators/integrals.hh>

#include "../problems/base.hh"
#include "base.hh"
#include "ipdg.hh"

namespace Dune {
namespace GDT {
namespace LinearElliptic {


/**
 * \todo add pattern() to StationaryContainerBasedDefaultDiscretization, avoid computation of local_pattern below
 */
template <class GridType,
          Backends spacebackend = Backends::gdt, // we only have local grid parts atm
          XT::LA::Backends la = XT::LA::default_sparse_backend,
          int pol = 1,
          class RangeFieldType = double,
          size_t dimRange = 1,
          LocalEllipticIpdgIntegrands::Method local_method = LocalEllipticIpdgIntegrands::default_method,
          LocalEllipticIpdgIntegrands::Method coupling_method = LocalEllipticIpdgIntegrands::default_method>
class BlockIpdgDiscretizer
{
  typedef IpdgDiscretizer<GridType,
                          XT::Grid::Layers::dd_subdomain,
                          spacebackend,
                          la,
                          pol,
                          RangeFieldType,
                          dimRange,
                          local_method>
      LocalDiscretizer;
  typedef typename LocalDiscretizer::DiscretizationType LocalDiscretizationType;
  typedef typename LocalDiscretizer::SpaceType LocalSpaceType;

public:
  typedef typename LocalDiscretizer::ProblemType ProblemType;
  typedef BlockSpace<LocalSpaceType> SpaceType;
  typedef typename LocalDiscretizer::MatrixType MatrixType;
  typedef typename LocalDiscretizer::VectorType VectorType;
  typedef StationaryContainerBasedDefaultDiscretization<ProblemType, SpaceType, MatrixType, VectorType, SpaceType>
      DiscretizationType;
  static const constexpr ChooseDiscretizer type = ChooseDiscretizer::block_ipdg;
  static const constexpr XT::LA::Backends la_backend = la;
  static const int polOrder = pol;

  static std::string static_id()
  {
    return std::string("gdt.linearelliptic.discretization.block-ipdg.order_")
           + Dune::XT::Common::to_string(int(polOrder)); // int() needed, otherwise we get a linker error
  }

  static DiscretizationType
  discretize(XT::Grid::GridProvider<GridType, XT::Grid::DD::SubdomainGrid<GridType>>& grid_provider,
             const ProblemType& problem,
             const int /*level*/ = 0)
  {
    auto logger = XT::Common::TimedLogger().get(static_id());
    const auto& dd_grid = grid_provider.dd_grid();
    logger.info() << "Creating " << dd_grid.size() << " local discretizations... " << std::endl;
    const LinearElliptic::ProblemBase<typename GridType::template Codim<0>::Entity,
                                      typename GridType::ctype,
                                      GridType::dimension,
                                      typename SpaceType::RangeFieldType,
                                      1>
        local_problem(problem.diffusion_factor(),
                      problem.diffusion_tensor(),
                      problem.force(),
                      problem.dirichlet(),
                      problem.neumann(),
                      XT::Common::Configuration(),
                      XT::Grid::allneumann_boundaryinfo_default_config());
    std::vector<LocalDiscretizationType> local_discretizations;
    for (size_t ss = 0; ss < dd_grid.size(); ++ss)
      local_discretizations.emplace_back(
          LocalDiscretizer::discretize(grid_provider, local_problem, XT::Common::numeric_cast<int>(ss)));

    logger.info() << "Creating space... " << std::endl;
    std::vector<std::shared_ptr<const LocalSpaceType>> local_spaces(dd_grid.size());
    for (size_t ss = 0; ss < dd_grid.size(); ++ss)
      local_spaces[ss] = std::make_shared<const LocalSpaceType>(local_discretizations[ss].test_space());
    SpaceType space(dd_grid, local_spaces);

    logger.info() << "Preparing container..." << std::endl;
    std::vector<XT::LA::SparsityPatternDefault> local_patterns(dd_grid.size());
    std::map<size_t, XT::LA::SparsityPatternDefault> boundary_patterns;
    std::vector<std::map<size_t, XT::LA::SparsityPatternDefault>> coupling_patterns(dd_grid.size());
    XT::LA::SparsityPatternDefault pattern(space.mapper().size());
    for (size_t ss = 0; ss < dd_grid.size(); ++ss) {
      local_patterns[ss] = local_spaces[ss]->compute_face_and_volume_pattern();
      copy_local_to_global(local_patterns[ss], space, ss, ss, pattern);
      if (dd_grid.boundary(ss)) {
        boundary_patterns[ss] = local_spaces[ss]->compute_volume_pattern(dd_grid.boundary_grid_view(ss));
        copy_local_to_global(boundary_patterns[ss], space, ss, ss, pattern);
      }
      coupling_patterns[ss][ss] =
          local_spaces[ss]->compute_volume_pattern(); // this is too much, but we are relazy here
      for (auto&& nn : dd_grid.neighborsOf(ss)) {
        if (ss < nn) { //  each coupling only once
          coupling_patterns[ss][nn] =
              local_spaces[ss]->compute_face_pattern(dd_grid.coupling_grid_view(ss, nn), *local_spaces[nn]);
          coupling_patterns[nn][ss] =
              local_spaces[nn]->compute_face_pattern(dd_grid.coupling_grid_view(nn, ss), *local_spaces[ss]);
          copy_local_to_global(coupling_patterns[ss][ss], space, ss, ss, pattern);
          copy_local_to_global(coupling_patterns[ss][nn], space, ss, nn, pattern);
          copy_local_to_global(coupling_patterns[nn][ss], space, nn, ss, pattern);
        }
      }
    }
    pattern.sort();
    MatrixType system_matrix(space.mapper().size(), space.mapper().size(), pattern);
    VectorType rhs_vector(space.mapper().size());
    for (size_t ss = 0; ss < dd_grid.size(); ++ss) {
      copy_local_to_global(local_discretizations[ss].system_matrix(), local_patterns[ss], space, ss, ss, system_matrix);
      copy_local_to_global(local_discretizations[ss].rhs_vector(), space, ss, rhs_vector);
    }

    logger.info() << "Assembling coupling contributions..." << std::endl;
    for (size_t ss = 0; ss < dd_grid.size(); ++ss) {
      for (auto&& nn : dd_grid.neighborsOf(ss)) {
        if (ss < nn) { //  each coupling only once
          MatrixType inside_inside_matrix(
              local_spaces[ss]->mapper().size(), local_spaces[ss]->mapper().size(), coupling_patterns[ss][ss]);
          MatrixType outside_outside_matrix(
              local_spaces[nn]->mapper().size(), local_spaces[nn]->mapper().size(), coupling_patterns[nn][nn]);
          MatrixType inside_outside_matrix(
              local_spaces[ss]->mapper().size(), local_spaces[nn]->mapper().size(), coupling_patterns[ss][nn]);
          MatrixType outside_inside_matrix(
              local_spaces[nn]->mapper().size(), local_spaces[ss]->mapper().size(), coupling_patterns[nn][ss]);
          auto coupling_grid_part = dd_grid.coupling_grid_view(ss, nn);
          // put all of this into a coupling operator
          SystemAssembler<LocalSpaceType, decltype(coupling_grid_part)> coupling_assembler(
              coupling_grid_part, *local_spaces[ss], *local_spaces[ss], *local_spaces[nn], *local_spaces[nn]);
          typedef XT::Grid::extract_intersection_t<decltype(coupling_grid_part)> IntersectionType;
          typedef LocalEllipticIpdgIntegrands::Inner<typename ProblemType::DiffusionFactorType,
                                                     typename ProblemType::DiffusionTensorType,
                                                     coupling_method>
              CouplingIntegrandType;
          LocalCouplingIntegralOperator<CouplingIntegrandType,
                                        typename LocalSpaceType::BaseFunctionSetType,
                                        IntersectionType>
              local_coupling_operator(problem.diffusion_factor(), problem.diffusion_tensor());
          coupling_assembler.append(local_coupling_operator,
                                    inside_inside_matrix,
                                    outside_outside_matrix,
                                    inside_outside_matrix,
                                    outside_inside_matrix);
          coupling_assembler.assemble();
          copy_local_to_global(inside_inside_matrix, coupling_patterns[ss][ss], space, ss, ss, system_matrix);
          copy_local_to_global(outside_outside_matrix, coupling_patterns[nn][nn], space, nn, nn, system_matrix);
          copy_local_to_global(inside_outside_matrix, coupling_patterns[ss][nn], space, ss, nn, system_matrix);
          copy_local_to_global(outside_inside_matrix, coupling_patterns[nn][ss], space, nn, ss, system_matrix);
        }
      }
    }

    logger.info() << "Assembling boundary contributions..." << std::endl;
    for (size_t ss = 0; ss < dd_grid.size(); ++ss) {
      if (dd_grid.boundary(ss)) {
        auto boundary_grid_part = dd_grid.boundary_grid_view(ss);
        MatrixType boundary_matrix(
            local_spaces[ss]->mapper().size(), local_spaces[ss]->mapper().size(), boundary_patterns[ss]);
        VectorType boundary_vector(local_spaces[ss]->mapper().size());
        SystemAssembler<LocalSpaceType, decltype(boundary_grid_part)> boundary_assembler(*local_spaces[ss],
                                                                                         boundary_grid_part);
        typedef XT::Grid::extract_intersection_t<decltype(boundary_grid_part)> IntersectionType;
        auto boundary_info = XT::Grid::BoundaryInfoFactory<IntersectionType>::create(problem.boundary_info_cfg());
        typedef LocalEllipticIpdgIntegrands::BoundaryLHS<typename ProblemType::DiffusionFactorType,
                                                         typename ProblemType::DiffusionTensorType,
                                                         coupling_method>
            BoundaryLhsIntegrandType;
        LocalBoundaryIntegralOperator<BoundaryLhsIntegrandType,
                                      typename LocalSpaceType::BaseFunctionSetType,
                                      IntersectionType>
            local_boundary_operator(problem.diffusion_factor(), problem.diffusion_tensor());
        boundary_assembler.append(
            local_boundary_operator,
            boundary_matrix,
            new XT::Grid::ApplyOn::DirichletIntersections<decltype(boundary_grid_part)>(*boundary_info));
        typedef LocalEllipticIpdgIntegrands::BoundaryRHS<typename ProblemType::FunctionType,
                                                         typename ProblemType::DiffusionFactorType,
                                                         typename ProblemType::DiffusionTensorType,
                                                         coupling_method>
            BoundaryRhsIntegrandType;
        LocalFaceIntegralFunctional<BoundaryRhsIntegrandType,
                                    typename LocalSpaceType::BaseFunctionSetType,
                                    IntersectionType>
            local_boundary_functional(problem.dirichlet(), problem.diffusion_factor(), problem.diffusion_tensor());
        LocalFaceFunctionalAssembler<LocalSpaceType, IntersectionType, VectorType> local_boundary_functional_assembler(
            local_boundary_functional);
        boundary_assembler.append(
            local_boundary_functional_assembler,
            boundary_vector,
            new XT::Grid::ApplyOn::DirichletIntersections<decltype(boundary_grid_part)>(*boundary_info));
        boundary_assembler.assemble();
        copy_local_to_global(boundary_matrix, boundary_patterns[ss], space, ss, ss, system_matrix);
        copy_local_to_global(boundary_vector, space, ss, rhs_vector);
      }
    }

    // create the discretization (no copy of the containers done here, bc. of cow)
    return DiscretizationType(problem, space, system_matrix, rhs_vector);
  } // ... discretize(...)

private:
  static void copy_local_to_global(const XT::LA::SparsityPatternDefault& local_pattern,
                                   const SpaceType& space,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   XT::LA::SparsityPatternDefault& global_pattern)
  {
    const auto& mapper = space.mapper();
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = mapper.mapToGlobal(test_subdomain, local_ii);
      const auto& local_rows = local_pattern.inner(local_ii);
      for (const auto& local_jj : local_rows) {
        const size_t global_jj = mapper.mapToGlobal(ansatz_subdomain, local_jj);
        global_pattern.insert(global_ii, global_jj);
      }
    }
  } // ... copy_local_to_global(...)

  static void copy_local_to_global(const MatrixType& local_matrix,
                                   const XT::LA::SparsityPatternDefault& local_pattern,
                                   const SpaceType& space,
                                   const size_t test_subdomain,
                                   const size_t ansatz_subdomain,
                                   MatrixType& global_matrix)
  {
    for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
      const size_t global_ii = space.mapper().mapToGlobal(test_subdomain, local_ii);
      for (const size_t& local_jj : local_pattern.inner(local_ii)) {
        const size_t global_jj = space.mapper().mapToGlobal(ansatz_subdomain, local_jj);
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
      }
    }
  } // ... copy_local_to_global_matrix(...)

  static void copy_local_to_global(const VectorType& local_vector,
                                   const SpaceType& space,
                                   const size_t subdomain,
                                   VectorType& global_vector)
  {
    for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
      const size_t global_ii = space.mapper().mapToGlobal(subdomain, local_ii);
      global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
    }
  } // ... copy_local_to_global_matrix(...)
}; // class BlockIpdgDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BLOCK_IPDG_HH
