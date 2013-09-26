// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH
#define DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localevaluation/sipdg.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretefunction/difference.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>

namespace Example {


template< class GridPartType >
class DiscretizationBase
{
public:
  static const unsigned int             dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int dimRange = 1;
  typedef double            RangeFieldType;

  typedef Dune::Stuff::FunctionInterface< DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::IntersectionType > BoundaryInfoType;

  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix< RangeFieldType >  MatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector< RangeFieldType >           VectorType;

  DiscretizationBase(const GridPartType& gp,
                     const BoundaryInfoType& info,
                     const FunctionType& diff,
                     const FunctionType& forc,
                     const FunctionType& dir,
                     const FunctionType& neu)
    : grid_part_(gp)
    , boundary_info_(info)
    , diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neu)
  {}

  const GridPartType& grid_part() const
  {
    return grid_part_;
  }

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  const FunctionType& diffusion() const
  {
    return diffusion_;
  }

  const FunctionType& force() const
  {
    return force_;
  }

  const FunctionType& dirichlet() const
  {
    return dirichlet_;
  }

  const FunctionType& neumann() const
  {
    return neumann_;
  }

  template< class ReferenceSolutionType, class DiscreteSolutionType >
  std::vector< RangeFieldType > compute_errors(const ReferenceSolutionType& referencec_solution,
                                               const DiscreteSolutionType& discrete_solution) const
  {
    using namespace Dune;
    using namespace Dune::GDT;
    // checks
    static_assert(std::is_base_of< Stuff::LocalizableFunction, ReferenceSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    static_assert(std::is_base_of< Stuff::LocalizableFunction, DiscreteSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    // difference
    typedef DiscreteFunction::Difference< ReferenceSolutionType, DiscreteSolutionType > DifferenceType;
    const DifferenceType difference(referencec_solution, discrete_solution);
    // L2 error
    ProductOperator::L2< GridPartType > l2_product_operator(grid_part_);
    const RangeFieldType l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    // H1 error
    ProductOperator::H1< GridPartType > h1_product_operator(grid_part_);
    const RangeFieldType h1_error = std::sqrt(h1_product_operator.apply2(difference, difference));
    return { l2_error, h1_error };
  }

private:
  const GridPartType& grid_part_;
  const BoundaryInfoType& boundary_info_;
  const FunctionType& diffusion_;
  const FunctionType& force_;
  const FunctionType& dirichlet_;
  const FunctionType& neumann_;
}; // class DiscretizationBase


template< class GridPartType, int polOrder >
class CGDiscretization
  : public DiscretizationBase< GridPartType >
{
  typedef DiscretizationBase< GridPartType > BaseType;
public:
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  typedef typename BaseType::FunctionType   FunctionType;

  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;

  typedef Dune::GDT::ContinuousLagrangeSpace::FemWrapper< GridPartType, polOrder, RangeFieldType, dimRange >  SpaceType;
  typedef Dune::GDT::DiscreteFunctionDefault< SpaceType, VectorType > DiscreteFunctionType;

  CGDiscretization(const GridPartType& gp,
                   const BoundaryInfoType& info,
                   const FunctionType& diff,
                   const FunctionType& forc,
                   const FunctionType& dir,
                   const FunctionType& neu)
    : BaseType(gp, info, diff, forc, dir, neu)
    , space_(BaseType::grid_part())
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const std::string id() const
  {
    return "cg." + Dune::Stuff::Common::toString(polOrder);
  }

  std::shared_ptr< DiscreteFunctionType > solve() const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    // left hand side
    // * elliptic diffusion operator
    typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< FunctionType > >  EllipticOperatorType;
    const EllipticOperatorType diffusion_operator(BaseType::diffusion());
    // * right hand side
    //   * L2 force functional
    typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > L2VolumeFunctionalType;
    const L2VolumeFunctionalType force_functional(BaseType::force());
    //   * L2 neumann functional
    typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > L2FaceFunctionalType;
    const L2FaceFunctionalType neumann_functional(BaseType::neumann());

    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType force_vector(space_.mapper().size());
    auto dirichlet_vector = std::make_shared< VectorType >(space_.mapper().size());
    VectorType neumann_vector(space_.mapper().size());
    VectorType rhs_vector(space_.mapper().size());
    auto solution = std::make_shared< DiscreteFunctionType >(space_,
                                                             std::make_shared< VectorType >(space_.mapper().size()));

    // * dirichlet boundary values
    DiscreteFunctionType dirichlet_projection(space_, dirichlet_vector, "dirichlet");

    typedef ProjectionOperator::Dirichlet< GridPartType > DirichletProjectionOperatorType;
    const DirichletProjectionOperatorType dirichlet_projection_operator(BaseType::grid_part(),
                                                                        BaseType::boundary_info());
    dirichlet_projection_operator.apply(BaseType::dirichlet(), dirichlet_projection);

    // * local matrix assembler
    typedef LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalEllipticOperatorMatrixAssemblerType;
    const LocalEllipticOperatorMatrixAssemblerType diffusion_matrix_assembler(diffusion_operator);
    // * local vector assemblers
    //   * force vector
    typedef LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalL2VolumeFunctionalVectorAssemblerType;
    const LocalL2VolumeFunctionalVectorAssemblerType force_vector_assembler(force_functional);
  //   * neumann vector
  typedef LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalL2FaceFunctionalVectorAssemblerType;
  const LocalL2FaceFunctionalVectorAssemblerType neumann_vector_assembler(neumann_functional);
    // * system assembler
    typedef SystemAssembler< SpaceType > SystemAssemblerType;
    SystemAssemblerType system_assembler(space_);
    system_assembler.addLocalAssembler(diffusion_matrix_assembler, system_matrix);
    system_assembler.addLocalAssembler(force_vector_assembler, force_vector);
    system_assembler.addLocalAssembler(neumann_vector_assembler,
                                      typename SystemAssemblerType::AssembleOnNeumann(BaseType::boundary_info()),
                                      neumann_vector);
    system_assembler.assemble();

    Constraints::Dirichlet< typename GridPartType::IntersectionType,
        RangeFieldType > dirichlet_constraints(BaseType::boundary_info(),
                                                                   space_.mapper().maxNumDofs(),
                                                                   space_.mapper().maxNumDofs());
    rhs_vector.backend() = force_vector.backend()
        + neumann_vector.backend()
        - system_matrix.backend()*dirichlet_vector->backend();
    system_assembler.addLocalConstraints(dirichlet_constraints, system_matrix);
    system_assembler.addLocalConstraints(dirichlet_constraints, rhs_vector);
    system_assembler.applyConstraints();

    typedef typename Dune::Stuff::LA::BicgstabILUTSolver< MatrixType, VectorType > LinearSolverType;
    auto linear_solver_settings = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-16";
    LinearSolverType linear_solver;
    const size_t failure = linear_solver.apply(system_matrix,
                                               rhs_vector,
                                               *(solution->vector()),
                                               linear_solver_settings);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver reported a problem!");
    if (solution->vector()->size() != space_.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver produced a solution of wrong size (is "
                 << solution->vector()->size() << ", should be " << space_.mapper().size() << ")!");
    solution->vector()->backend() += dirichlet_vector->backend();
    return solution;
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename) const
  {
    space_.visualize(vector, filename);
  }

private:
  const SpaceType space_;
}; // class CGDiscretization


template< class GridPartType, int polOrder >
class SIPDGDiscretization
  : public DiscretizationBase< GridPartType >
{
  typedef DiscretizationBase< GridPartType > BaseType;
public:
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  typedef typename BaseType::FunctionType   FunctionType;

  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;

  typedef Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<  GridPartType,
                                                                            polOrder,
                                                                            RangeFieldType,
                                                                            dimRange > SpaceType;
  typedef Dune::GDT::DiscreteFunctionDefault< SpaceType, VectorType > DiscreteFunctionType;

  SIPDGDiscretization(const GridPartType& gp,
                      const BoundaryInfoType& info,
                      const FunctionType& diff,
                      const FunctionType& forc,
                      const FunctionType& dir,
                      const FunctionType& neu)
    : BaseType(gp, info, diff, forc, dir, neu)
    , space_(BaseType::grid_part())
    , penalty_factor_(1.0)
    , integration_order_(4)
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const std::string id() const
  {
    return "sipdg." + Dune::Stuff::Common::toString(polOrder);
  }

  std::shared_ptr< DiscreteFunctionType > solve() const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    // container
    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType rhs_vector(space_.mapper().size());
    auto solution = std::make_shared< DiscreteFunctionType >(space_,
                                                             std::make_shared< VectorType >(space_.mapper().size()));

    // walk the grid
    const auto entity_it_end = BaseType::grid_part().template end< 0 >();
    for (auto entity_it = BaseType::grid_part().template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      // entity
      const auto& entity = *entity_it;
      const auto test_basis_entity = space_.baseFunctionSet(entity);
      const auto ansatz_basis_entity = space_.baseFunctionSet(entity);
      const auto local_diffusion_entity = BaseType::diffusion().localFunction(entity);
      const auto local_force = BaseType::force().localFunction(entity);
      const auto local_dirichlet = BaseType::dirichlet().localFunction(entity);
      const auto local_neumann = BaseType::neumann().localFunction(entity);
      const auto global_test_indices_entity = space_.mapper().globalIndices(entity);
      const auto global_ansatz_indices_entity = space_.mapper().globalIndices(entity);

      // prepare local storage
      DynamicMatrix< RangeFieldType > lhs_volume_integrals(test_basis_entity.size(),ansatz_basis_entity.size(), 0);
      DynamicVector< RangeFieldType > rhs_volume_integrals(test_basis_entity.size(), 0);
      // do a volume quadrature
      const auto& volume_quadratue = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                         2*integration_order_ + 1);
      for (const auto& quadrature_point : volume_quadratue) {
        // evaluate
        const auto quadrature_weight = quadrature_point.weight();
        const auto local_point_entity = quadrature_point.position();
        const auto integration_element = entity.geometry().integrationElement(local_point_entity);
        const auto diffusion_value = local_diffusion_entity.evaluate(local_point_entity);
        const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
        const auto test_base_gradients = test_basis_entity.jacobian(local_point_entity);
        const auto ansatz_base_gradients = ansatz_basis_entity.jacobian(local_point_entity);
        const auto force_value = local_force.evaluate(local_point_entity);
        // compute integrals
        // * lsh
        for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
          for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
            lhs_volume_integrals[ii][jj] += quadrature_weight * integration_element
                                            * diffusion_value * (ansatz_base_gradients[jj][0] * test_base_gradients[ii][0]);
          }
        }
        // * rhs
        for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
          rhs_volume_integrals[ii] += quadrature_weight * integration_element * (force_value * test_base_values[ii]);
        }
      } // do a volume quadrature

      // map
      // * lhs
      for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
        for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
          const auto& global_ii = global_test_indices_entity[ii];
          const auto& global_jj = global_ansatz_indices_entity[jj];
          system_matrix.add(global_ii, global_jj, lhs_volume_integrals[ii][jj]);
        }
      }
      // * rhs
      for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
        const auto& global_ii = global_test_indices_entity[ii];
        rhs_vector.add(global_ii, rhs_volume_integrals[ii]);
      }

      // walk the intersections
      const auto intersection_it_end = BaseType::grid_part().iend(entity);
      for (auto intersection_it = BaseType::grid_part().ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          // we are on an inner intersection
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // work only once on each intersection
          if (BaseType::grid_part().indexSet().index(entity) < BaseType::grid_part().indexSet().index(neighbour)) {
            // neighbour
            const auto test_basis_neighbour = space_.baseFunctionSet(neighbour);
            const auto ansatz_basis_neighbour = space_.baseFunctionSet(neighbour);
            const auto local_diffusion_neighbour = BaseType::diffusion().localFunction(neighbour);
            const auto global_test_indices_neighbour = space_.mapper().globalIndices(neighbour);
            const auto global_ansatz_indices_neighbour = space_.mapper().globalIndices(neighbour);
            // compute local face weight
            RangeFieldType local_face_weight = /*0;
            {
              assert(dimDomain > 1);
              assert(intersection.geometry().corners() == 2);
              const auto corner_difference = intersection.geometry().corner(0) - intersection.geometry().corner(1);
              local_face_weight = corner_difference.two_norm();
            }*/ intersection.geometry().volume();
            // prepare local storage
            DynamicMatrix< RangeFieldType > integrals_entity_entity(test_basis_entity.size(),
                                                                    test_basis_entity.size(),
                                                                    0);
            DynamicMatrix< RangeFieldType > integrals_entity_neighbour(test_basis_entity.size(),
                                                                       test_basis_neighbour.size(),
                                                                       0);
            DynamicMatrix< RangeFieldType > integrals_neighbour_entity(test_basis_neighbour.size(),
                                                                       test_basis_entity.size(),
                                                                       0);
            DynamicMatrix< RangeFieldType > integrals_neighbour_neighbour(test_basis_neighbour.size(),
                                                                          test_basis_neighbour.size(),
                                                                          0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto unit_outer_normal = intersection.unitOuterNormal(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto local_point_neighbour = intersection.geometryInOutside().global(local_point_intersection);
              const auto diffusion_value_entity = local_diffusion_entity.evaluate(local_point_entity);
              const auto test_base_values_entity = test_basis_entity.evaluate(local_point_entity);
              const auto test_base_gradients_entity = test_basis_entity.jacobian(local_point_entity);
              const auto ansatz_base_values_entity = ansatz_basis_entity.evaluate(local_point_entity);
              const auto ansatz_base_gradients_entity = ansatz_basis_entity.jacobian(local_point_entity);
              const auto diffusion_value_neighbour = local_diffusion_neighbour.evaluate(local_point_neighbour);
              const auto test_base_values_neighbour = test_basis_neighbour.evaluate(local_point_neighbour);
              const auto test_base_gradients_neighbour = test_basis_neighbour.jacobian(local_point_neighbour);
              const auto ansatz_base_values_neighbour = ansatz_basis_neighbour.evaluate(local_point_neighbour);
              const auto ansatz_base_gradients_neighbour = ansatz_basis_neighbour.jacobian(local_point_neighbour);
              // compute integrals
              // * entity / entity
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (-0.5 * diffusion_value_entity * (ansatz_base_gradients_entity[jj][0] * unit_outer_normal) * test_base_values_entity[ii]);
                  // symmetry term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (-0.5 * ansatz_base_values_entity[jj] * diffusion_value_entity * (test_base_gradients_entity[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (((8.0 * penalty_factor_) / local_face_weight) * ansatz_base_values_entity[jj] * test_base_values_entity[ii]);
                }
              }
              // * entity / neighbour
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                  // consistency term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (-0.5 * diffusion_value_neighbour * (ansatz_base_gradients_neighbour[jj][0] * unit_outer_normal) * test_base_values_entity[ii]);
                  // symmetry term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (0.5 * ansatz_base_values_neighbour[jj] * diffusion_value_entity * (test_base_gradients_entity[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * ((8.0 * penalty_factor_) / local_face_weight) * ansatz_base_values_neighbour[jj] * test_base_values_entity[ii]);
                }
              }
              // * neighbour / entity
              for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (0.5 * diffusion_value_entity * (ansatz_base_gradients_entity[jj][0] * unit_outer_normal) * test_base_values_neighbour[ii]);
                  // symmetry term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (-0.5 * ansatz_base_values_entity[jj] * diffusion_value_neighbour * (test_base_gradients_neighbour[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * ((8.0 * penalty_factor_) / local_face_weight) * ansatz_base_values_entity[jj] * test_base_values_neighbour[ii]);
                }
              }
              // * neighbour / neighbour
              for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                  // consistency term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (0.5 * diffusion_value_neighbour * (ansatz_base_gradients_neighbour[jj][0] * unit_outer_normal) * test_base_values_neighbour[ii]);
                  // symmetry term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (0.5 * ansatz_base_values_neighbour[jj] * diffusion_value_neighbour * (test_base_gradients_neighbour[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (((8.0 * penalty_factor_) / local_face_weight) * ansatz_base_values_neighbour[jj] * test_base_values_neighbour[ii]);
                }
              }
            } // do a face quadrature

            // map
            // * entity / entity
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, integrals_entity_entity[ii][jj]);
              }
            }
            // * entity / neighbour
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_neighbour[jj];
                system_matrix.add(global_ii, global_jj, integrals_entity_neighbour[ii][jj]);
              }
            }
            // * neighbour / entity
            for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_neighbour[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, integrals_neighbour_entity[ii][jj]);
              }
            }
            // * neighbour / neighbour
            for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                const auto& global_ii = global_test_indices_neighbour[ii];
                const auto& global_jj = global_ansatz_indices_neighbour[jj];
                system_matrix.add(global_ii, global_jj, integrals_neighbour_neighbour[ii][jj]);
              }
            }
          } // work only once on each intersection
        } else if (intersection.boundary() && !intersection.neighbor()) {
          // we are on a boundary intersection
          // compute local face weight
          RangeFieldType local_face_weight = /*0;
          {
            assert(dimDomain > 1);
            assert(intersection.geometry().corners() == 2);
            const auto corner_difference = intersection.geometry().corner(0) - intersection.geometry().corner(1);
            local_face_weight = corner_difference.two_norm();
          }*/ intersection.geometry().volume();

          // handle dirichlet
          if (BaseType::boundary_info().dirichlet(intersection)) {
            // prepare local storage
            DynamicMatrix< RangeFieldType > lhs_integrals(test_basis_entity.size(), test_basis_entity.size(), 0);
            DynamicVector< RangeFieldType > rhs_face_integrals(test_basis_entity.size(), 0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto unit_outer_normal = intersection.unitOuterNormal(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto diffusion_value = local_diffusion_entity.evaluate(local_point_entity);
              const auto dirichlet_value = local_dirichlet.evaluate(local_point_entity);
              const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
              const auto test_base_gradients = test_basis_entity.jacobian(local_point_entity);
              const auto ansatz_base_values = ansatz_basis_entity.evaluate(local_point_entity);
              const auto ansatz_base_gradients = ansatz_basis_entity.jacobian(local_point_entity);
              // compute integrals
              // * lhs
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * diffusion_value * (ansatz_base_gradients[jj][0] * unit_outer_normal) * test_base_values[ii]);
                  // symmetry term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * ansatz_base_values[jj] * diffusion_value * (test_base_gradients[ii][0] * unit_outer_normal));
                  // penalty term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (((14.0 * penalty_factor_) / local_face_weight) * ansatz_base_values[jj] * test_base_values[ii]);
                }
              }
              // * rhs
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                // dirichlet symmetry term
                rhs_face_integrals[ii] -= quadrature_weight * integration_element
                    * (dirichlet_value * diffusion_value * (test_base_gradients[ii][0] * unit_outer_normal));
                // dirichlet penalty term
                rhs_face_integrals[ii] += quadrature_weight * integration_element
                    * (((14.0 * penalty_factor_) / local_face_weight) * dirichlet_value * test_base_values[ii]);
              }
            } // do a face quadrature

            // map
            // * lhs
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, lhs_integrals[ii][jj]);
              }
            }
            // * rhs
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              const auto& global_ii = global_test_indices_entity[ii];
              rhs_vector.add(global_ii, rhs_face_integrals[ii]);
            }
          } else if (BaseType::boundary_info().neumann(intersection)) {
            // handle neumann
            // prepare local storage
            DynamicVector< RangeFieldType > rhs_face_integrals(test_basis_entity.size(), 0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto neumann_value = local_neumann.evaluate(local_point_entity);
              const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
              // compute integrals
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii)
                rhs_face_integrals[ii] -= quadrature_weight * integration_element
                    * (neumann_value * test_base_values[ii]);
            } // do a face quadrature

            // map
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              const auto& global_ii = global_test_indices_entity[ii];
              rhs_vector.add(global_ii, rhs_face_integrals[ii]);
            }
          } else
            DUNE_THROW(InvalidStateException, "Unknown type of boundary information encountered!");
        } else
          DUNE_THROW(InvalidStateException, "Unknown type of intersection encountered!");
      } // walk the intersections
    } // walk the grid

    // solve
    typedef typename Dune::Stuff::LA::BicgstabILUTSolver< MatrixType, VectorType > LinearSolverType;
    auto linear_solver_settings = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-16";
    LinearSolverType linear_solver;
    const size_t failure = linear_solver.apply(system_matrix,
                                               rhs_vector,
                                               *(solution->vector()),
                                               linear_solver_settings);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver reported a problem!");
    if (solution->vector()->size() != space_.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver produced a solution of wrong size (is "
                 << solution->vector()->size() << ", should be " << space_.mapper().size() << ")!");

    return solution;
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename) const
  {
    space_.visualize(vector, filename);
  }

private:
  const SpaceType space_;
  const RangeFieldType penalty_factor_;
  const size_t integration_order_;
}; // class SIPDGDiscretization


template< class GridPartType, int polOrder >
class NewSIPDGDiscretization
  : public DiscretizationBase< GridPartType >
{
  typedef DiscretizationBase< GridPartType > BaseType;
public:
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  typedef typename BaseType::FunctionType   FunctionType;

  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;

  typedef Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<  GridPartType,
                                                                            polOrder,
                                                                            RangeFieldType,
                                                                            dimRange > SpaceType;
  typedef Dune::GDT::DiscreteFunctionDefault< SpaceType, VectorType > DiscreteFunctionType;

  NewSIPDGDiscretization(const GridPartType& gp,
                         const BoundaryInfoType& info,
                         const FunctionType& diff,
                         const FunctionType& forc,
                         const FunctionType& dir,
                         const FunctionType& neu)
    : BaseType(gp, info, diff, forc, dir, neu)
    , space_(BaseType::grid_part())
    , beta_(1.0)
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const std::string id() const
  {
    return "newsipdg." + Dune::Stuff::Common::toString(polOrder);
  }

  std::shared_ptr< DiscreteFunctionType > solve() const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    // container
    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType rhs_vector(space_.mapper().size());
    auto solution = std::make_shared< DiscreteFunctionType >(space_,
                                                             std::make_shared< VectorType >(space_.mapper().size()));

    typedef SystemAssembler< SpaceType > SystemAssemblerType;
    SystemAssemblerType systemAssembler(space_);

    // volume terms
    // * lhs
    typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< FunctionType > > EllipticOperatorType;
    const EllipticOperatorType                                  ellipticOperator(BaseType::diffusion());
    const LocalAssembler::Codim0Matrix< EllipticOperatorType >  diffusionMatrixAssembler(ellipticOperator);
    systemAssembler.addLocalAssembler(diffusionMatrixAssembler, system_matrix);
    // * rhs
    typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > ForceFunctionalType;
    const ForceFunctionalType                                 forceFunctional(BaseType::force());
    const LocalAssembler::Codim0Vector< ForceFunctionalType > forceVectorAssembler(forceFunctional);
    systemAssembler.addLocalAssembler(forceVectorAssembler, rhs_vector);
    // inner face terms
    typedef LocalOperator::Codim1CouplingIntegral< LocalEvaluation::SIPDG::Inner< FunctionType > > CouplingOperatorType;
    const CouplingOperatorType                                          couplingOperator(BaseType::diffusion(), beta_);
    const LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType >  couplingMatrixAssembler(couplingOperator);
    systemAssembler.addLocalAssembler(couplingMatrixAssembler,
                                      typename SystemAssemblerType::AssembleOnInnerPrimally(),
                                      system_matrix);
    // dirichlet boundary face terms
    // * lhs
    typedef LocalOperator::Codim1BoundaryIntegral< LocalEvaluation::SIPDG::BoundaryLHS< FunctionType > >
        DirichletOperatorType;
    const DirichletOperatorType                                         dirichletOperator(BaseType::diffusion(), beta_);
    const LocalAssembler::Codim1BoundaryMatrix< DirichletOperatorType > dirichletMatrixAssembler(dirichletOperator);
    systemAssembler.addLocalAssembler(dirichletMatrixAssembler,
                                      typename SystemAssemblerType::AssembleOnDirichlet(BaseType::boundary_info()),
                                      system_matrix);
    // * rhs
    typedef LocalFunctional::Codim1Integral< LocalEvaluation::SIPDG::BoundaryRHS< FunctionType, FunctionType > >
        DirichletFunctionalType;
    const DirichletFunctionalType                                 dirichletFunctional(BaseType::diffusion(),
                                                                                      BaseType::dirichlet(),
                                                                                      beta_);
    const LocalAssembler::Codim1Vector< DirichletFunctionalType > dirichletVectorAssembler(dirichletFunctional);
    systemAssembler.addLocalAssembler(dirichletVectorAssembler,
                                      typename SystemAssemblerType::AssembleOnDirichlet(BaseType::boundary_info()),
                                      rhs_vector);
    // neumann boundary face terms
    // * rhs
    typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > NeumannFunctionalType;
    const NeumannFunctionalType                                 neumannFunctional(BaseType::neumann());
    const LocalAssembler::Codim1Vector< NeumannFunctionalType > neumannVectorAssembler(neumannFunctional);
    systemAssembler.addLocalAssembler(neumannVectorAssembler,
                                      typename SystemAssemblerType::AssembleOnNeumann(BaseType::boundary_info()),
                                      rhs_vector);

    // do all the work
    systemAssembler.assemble();

    // solve
    typedef typename Dune::Stuff::LA::BicgstabILUTSolver< MatrixType, VectorType > LinearSolverType;
    auto linear_solver_settings = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-16";
    LinearSolverType linear_solver;
    const size_t failure = linear_solver.apply(system_matrix,
                                               rhs_vector,
                                               *(solution->vector()),
                                               linear_solver_settings);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver reported a problem!");
    if (solution->vector()->size() != space_.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver produced a solution of wrong size (is "
                 << solution->vector()->size() << ", should be " << space_.mapper().size() << ")!");

    return solution;
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename) const
  {
    space_.visualize(vector, filename);
  }

private:
  const SpaceType space_;
  const RangeFieldType beta_;
}; // class NewSIPDGDiscretization


template< class GridPartType, int polOrder >
class SWIPDGDiscretization
  : public DiscretizationBase< GridPartType >
{
  typedef DiscretizationBase< GridPartType > BaseType;
public:
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  typedef typename BaseType::FunctionType   FunctionType;

  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;

  typedef Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<  GridPartType,
                                                                            polOrder,
                                                                            RangeFieldType,
                                                                            dimRange > SpaceType;
  typedef Dune::GDT::DiscreteFunctionDefault< SpaceType, VectorType > DiscreteFunctionType;

  SWIPDGDiscretization(const GridPartType& gp,
                       const BoundaryInfoType& info,
                       const FunctionType& diff,
                       const FunctionType& forc,
                       const FunctionType& dir,
                       const FunctionType& neu)
    : BaseType(gp, info, diff, forc, dir, neu)
    , space_(BaseType::grid_part())
    , penalty_factor_(10.0)
    , integration_order_(4)
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const std::string id() const
  {
    return "swipdg." + Dune::Stuff::Common::toString(polOrder);
  }

  std::shared_ptr< DiscreteFunctionType > solve() const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    // container
    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType rhs_vector(space_.mapper().size());
    auto solution = std::make_shared< DiscreteFunctionType >(space_,
                                                             std::make_shared< VectorType >(space_.mapper().size()));

    // walk the grid
    const auto entity_it_end = BaseType::grid_part().template end< 0 >();
    for (auto entity_it = BaseType::grid_part().template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      // entity
      const auto& entity = *entity_it;
      const auto test_basis_entity = space_.baseFunctionSet(entity);
      const auto ansatz_basis_entity = space_.baseFunctionSet(entity);
      const auto local_diffusion_entity = BaseType::diffusion().localFunction(entity);
      const auto local_force = BaseType::force().localFunction(entity);
      const auto local_dirichlet = BaseType::dirichlet().localFunction(entity);
      const auto local_neumann = BaseType::neumann().localFunction(entity);
      const auto global_test_indices_entity = space_.mapper().globalIndices(entity);
      const auto global_ansatz_indices_entity = space_.mapper().globalIndices(entity);

      // prepare local storage
      DynamicMatrix< RangeFieldType > lhs_volume_integrals(test_basis_entity.size(),ansatz_basis_entity.size(), 0);
      DynamicVector< RangeFieldType > rhs_volume_integrals(test_basis_entity.size(), 0);
      // do a volume quadrature
      const auto& volume_quadratue = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                         2*integration_order_ + 1);
      for (const auto& quadrature_point : volume_quadratue) {
        // evaluate
        const auto quadrature_weight = quadrature_point.weight();
        const auto local_point_entity = quadrature_point.position();
        const auto integration_element = entity.geometry().integrationElement(local_point_entity);
        const auto diffusion_value = local_diffusion_entity.evaluate(local_point_entity);
        const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
        const auto test_base_gradients = test_basis_entity.jacobian(local_point_entity);
        const auto ansatz_base_gradients = ansatz_basis_entity.jacobian(local_point_entity);
        const auto force_value = local_force.evaluate(local_point_entity);
        // compute integrals
        // * lsh
        for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
          for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
            lhs_volume_integrals[ii][jj] += quadrature_weight * integration_element
                                            * diffusion_value * (ansatz_base_gradients[jj][0] * test_base_gradients[ii][0]);
          }
        }
        // * rhs
        for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
          rhs_volume_integrals[ii] += quadrature_weight * integration_element * (force_value * test_base_values[ii]);
        }
      } // do a volume quadrature

      // map
      // * lhs
      for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
        for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
          const auto& global_ii = global_test_indices_entity[ii];
          const auto& global_jj = global_ansatz_indices_entity[jj];
          system_matrix.add(global_ii, global_jj, lhs_volume_integrals[ii][jj]);
        }
      }
      // * rhs
      for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
        const auto& global_ii = global_test_indices_entity[ii];
        rhs_vector.add(global_ii, rhs_volume_integrals[ii]);
      }

      // walk the intersections
      const auto intersection_it_end = BaseType::grid_part().iend(entity);
      for (auto intersection_it = BaseType::grid_part().ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          // we are on an inner intersection
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // work only once on each intersection
          if (BaseType::grid_part().indexSet().index(entity) < BaseType::grid_part().indexSet().index(neighbour)) {
            // neighbour
            const auto test_basis_neighbour = space_.baseFunctionSet(neighbour);
            const auto ansatz_basis_neighbour = space_.baseFunctionSet(neighbour);
            const auto local_diffusion_neighbour = BaseType::diffusion().localFunction(neighbour);
            const auto global_test_indices_neighbour = space_.mapper().globalIndices(neighbour);
            const auto global_ansatz_indices_neighbour = space_.mapper().globalIndices(neighbour);
            // compute local face weight
            RangeFieldType local_face_weight = /*0;
            {
              assert(dimDomain > 1);
              assert(intersection.geometry().corners() == 2);
              const auto corner_difference = intersection.geometry().corner(0) - intersection.geometry().corner(1);
              local_face_weight = corner_difference.two_norm();
            }*/ intersection.geometry().volume();
            // prepare local storage
            DynamicMatrix< RangeFieldType > integrals_entity_entity(test_basis_entity.size(),
                                                                    test_basis_entity.size(),
                                                                    0);
            DynamicMatrix< RangeFieldType > integrals_entity_neighbour(test_basis_entity.size(),
                                                                       test_basis_neighbour.size(),
                                                                       0);
            DynamicMatrix< RangeFieldType > integrals_neighbour_entity(test_basis_neighbour.size(),
                                                                       test_basis_entity.size(),
                                                                       0);
            DynamicMatrix< RangeFieldType > integrals_neighbour_neighbour(test_basis_neighbour.size(),
                                                                          test_basis_neighbour.size(),
                                                                          0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto unit_outer_normal = intersection.unitOuterNormal(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto local_point_neighbour = intersection.geometryInOutside().global(local_point_intersection);
              const auto diffusion_value_entity = local_diffusion_entity.evaluate(local_point_entity);
              const auto test_base_values_entity = test_basis_entity.evaluate(local_point_entity);
              const auto test_base_gradients_entity = test_basis_entity.jacobian(local_point_entity);
              const auto ansatz_base_values_entity = ansatz_basis_entity.evaluate(local_point_entity);
              const auto ansatz_base_gradients_entity = ansatz_basis_entity.jacobian(local_point_entity);
              const auto diffusion_value_neighbour = local_diffusion_neighbour.evaluate(local_point_neighbour);
              const auto test_base_values_neighbour = test_basis_neighbour.evaluate(local_point_neighbour);
              const auto test_base_gradients_neighbour = test_basis_neighbour.jacobian(local_point_neighbour);
              const auto ansatz_base_values_neighbour = ansatz_basis_neighbour.evaluate(local_point_neighbour);
              const auto ansatz_base_gradients_neighbour = ansatz_basis_neighbour.jacobian(local_point_neighbour);
              // compute weight and penalty
              const RangeFieldType weight_entity = diffusion_value_neighbour
                                                   / (diffusion_value_entity + diffusion_value_neighbour);
              const RangeFieldType weight_neighbour = diffusion_value_entity
                                                      / (diffusion_value_entity + diffusion_value_neighbour);
              const RangeFieldType penalty = (penalty_factor_ * 2.0 * diffusion_value_entity * diffusion_value_neighbour)
                                             / ((diffusion_value_entity + diffusion_value_neighbour) * local_face_weight);
              // compute integrals
              // * entity / entity
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (-weight_entity * diffusion_value_entity * (ansatz_base_gradients_entity[jj][0] * unit_outer_normal) * test_base_values_entity[ii]);
                  // symmetry term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (-weight_entity * ansatz_base_values_entity[jj] * diffusion_value_entity * (test_base_gradients_entity[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_entity_entity[ii][jj] += quadrature_weight * integration_element
                      * (penalty * ansatz_base_values_entity[jj] * test_base_values_entity[ii]);
                }
              }
              // * entity / neighbour
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                  // consistency term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (-weight_neighbour * diffusion_value_neighbour * (ansatz_base_gradients_neighbour[jj][0] * unit_outer_normal) * test_base_values_entity[ii]);
                  // symmetry term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (weight_entity * ansatz_base_values_neighbour[jj] * diffusion_value_entity * (test_base_gradients_entity[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_entity_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * penalty * ansatz_base_values_neighbour[jj] * test_base_values_entity[ii]);
                }
              }
              // * neighbour / entity
              for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (weight_entity * diffusion_value_entity * (ansatz_base_gradients_entity[jj][0] * unit_outer_normal) * test_base_values_neighbour[ii]);
                  // symmetry term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (-weight_neighbour * ansatz_base_values_entity[jj] * diffusion_value_neighbour * (test_base_gradients_neighbour[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_neighbour_entity[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * penalty * ansatz_base_values_entity[jj] * test_base_values_neighbour[ii]);
                }
              }
              // * neighbour / neighbour
              for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                  // consistency term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (weight_neighbour * diffusion_value_neighbour * (ansatz_base_gradients_neighbour[jj][0] * unit_outer_normal) * test_base_values_neighbour[ii]);
                  // symmetry term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (weight_neighbour * ansatz_base_values_neighbour[jj] * diffusion_value_neighbour * (test_base_gradients_neighbour[ii][0] * unit_outer_normal));
                  // penalty term
                  integrals_neighbour_neighbour[ii][jj] += quadrature_weight * integration_element
                      * (penalty * ansatz_base_values_neighbour[jj] * test_base_values_neighbour[ii]);
                }
              }
            } // do a face quadrature

            // map
            // * entity / entity
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, integrals_entity_entity[ii][jj]);
              }
            }
            // * entity / neighbour
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_neighbour[jj];
                system_matrix.add(global_ii, global_jj, integrals_entity_neighbour[ii][jj]);
              }
            }
            // * neighbour / entity
            for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_neighbour[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, integrals_neighbour_entity[ii][jj]);
              }
            }
            // * neighbour / neighbour
            for (size_t ii = 0; ii < test_basis_neighbour.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_neighbour.size(); ++jj) {
                const auto& global_ii = global_test_indices_neighbour[ii];
                const auto& global_jj = global_ansatz_indices_neighbour[jj];
                system_matrix.add(global_ii, global_jj, integrals_neighbour_neighbour[ii][jj]);
              }
            }
          } // work only once on each intersection
        } else if (intersection.boundary() && !intersection.neighbor()) {
          // we are on a boundary intersection
          // compute local face weight
          RangeFieldType local_face_weight = /*0;
          {
            assert(dimDomain > 1);
            assert(intersection.geometry().corners() == 2);
            const auto corner_difference = intersection.geometry().corner(0) - intersection.geometry().corner(1);
            local_face_weight = corner_difference.two_norm();
          }*/ intersection.geometry().volume();

          // handle dirichlet
          if (BaseType::boundary_info().dirichlet(intersection)) {
            // prepare local storage
            DynamicMatrix< RangeFieldType > lhs_integrals(test_basis_entity.size(), test_basis_entity.size(), 0);
            DynamicVector< RangeFieldType > rhs_face_integrals(test_basis_entity.size(), 0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto unit_outer_normal = intersection.unitOuterNormal(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto diffusion_value = local_diffusion_entity.evaluate(local_point_entity);
              const auto dirichlet_value = local_dirichlet.evaluate(local_point_entity);
              const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
              const auto test_base_gradients = test_basis_entity.jacobian(local_point_entity);
              const auto ansatz_base_values = ansatz_basis_entity.evaluate(local_point_entity);
              const auto ansatz_base_gradients = ansatz_basis_entity.jacobian(local_point_entity);
              const RangeFieldType penalty = (penalty_factor_ * diffusion_value) / local_face_weight;
              // compute integrals
              // * lhs
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                  // consistency term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * diffusion_value * (ansatz_base_gradients[jj][0] * unit_outer_normal) * test_base_values[ii]);
                  // symmetry term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (-1.0 * ansatz_base_values[jj] * diffusion_value * (test_base_gradients[ii][0] * unit_outer_normal));
                  // penalty term
                  lhs_integrals[ii][jj] += quadrature_weight * integration_element
                      * (penalty * ansatz_base_values[jj] * test_base_values[ii]);
                }
              }
              // * rhs
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
                // dirichlet symmetry term
                rhs_face_integrals[ii] -= quadrature_weight * integration_element
                    * (dirichlet_value * diffusion_value * (test_base_gradients[ii][0] * unit_outer_normal));
                // dirichlet penalty term
                rhs_face_integrals[ii] += quadrature_weight * integration_element
                    * (penalty * dirichlet_value * test_base_values[ii]);
              }
            } // do a face quadrature

            // map
            // * lhs
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              for (size_t jj = 0; jj < ansatz_basis_entity.size(); ++jj) {
                const auto& global_ii = global_test_indices_entity[ii];
                const auto& global_jj = global_ansatz_indices_entity[jj];
                system_matrix.add(global_ii, global_jj, lhs_integrals[ii][jj]);
              }
            }
            // * rhs
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              const auto& global_ii = global_test_indices_entity[ii];
              rhs_vector.add(global_ii, rhs_face_integrals[ii]);
            }
          } else if (BaseType::boundary_info().neumann(intersection)) {
            // handle neumann
            // prepare local storage
            DynamicVector< RangeFieldType > rhs_face_integrals(test_basis_entity.size(), 0);
            // do a face quadrature
            const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                                  2*integration_order_ + 1);
            for (const auto& quadrature_point : face_quadrature) {
              // evaluate
              const auto quadrature_weight = quadrature_point.weight();
              const auto local_point_intersection = quadrature_point.position();
              const auto integration_element = intersection.geometry().integrationElement(local_point_intersection);
              const auto local_point_entity = intersection.geometryInInside().global(local_point_intersection);
              const auto neumann_value = local_neumann.evaluate(local_point_entity);
              const auto test_base_values = test_basis_entity.evaluate(local_point_entity);
              // compute integrals
              for (size_t ii = 0; ii < test_basis_entity.size(); ++ii)
                rhs_face_integrals[ii] -= quadrature_weight * integration_element
                    * (neumann_value * test_base_values[ii]);
            } // do a face quadrature

            // map
            for (size_t ii = 0; ii < test_basis_entity.size(); ++ii) {
              const auto& global_ii = global_test_indices_entity[ii];
              rhs_vector.add(global_ii, rhs_face_integrals[ii]);
            }
          } else
            DUNE_THROW(InvalidStateException, "Unknown type of boundary information encountered!");
        } else
          DUNE_THROW(InvalidStateException, "Unknown type of intersection encountered!");
      } // walk the intersections
    } // walk the grid

    // solve
    typedef typename Dune::Stuff::LA::BicgstabILUTSolver< MatrixType, VectorType > LinearSolverType;
    auto linear_solver_settings = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-16";
    LinearSolverType linear_solver;
    const size_t failure = linear_solver.apply(system_matrix,
                                               rhs_vector,
                                               *(solution->vector()),
                                               linear_solver_settings);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver reported a problem!");
    if (solution->vector()->size() != space_.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver produced a solution of wrong size (is "
                 << solution->vector()->size() << ", should be " << space_.mapper().size() << ")!");

    return solution;
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename) const
  {
    space_.visualize(vector, filename);
  }

private:
  const SpaceType space_;
  const RangeFieldType penalty_factor_;
  const size_t integration_order_;
}; // class SWIPDGDiscretization


} // namespace Example

#endif // DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH
