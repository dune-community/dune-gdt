// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_SWIPDG_HH
#define DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_SWIPDG_HH

#include "discretization-base.hh"

namespace Discretization {


template< class GridPartType, int polOrder >
class SymmetricWeightedInteriorPenaltyDG
  : public Base< GridPartType >
{
  typedef Base< GridPartType > BaseType;
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

  SymmetricWeightedInteriorPenaltyDG(const GridPartType& gp,
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
}; // class SymmetricWeightedInteriorPenaltyDG


} // namespace Discretization

#endif // DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_SWIPDG_HH
