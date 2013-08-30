// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH
#define DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH

#include <memory>
#include <type_traits>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretefunction/difference.hh>
#include <dune/gdt/operator/projections.hh>

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

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::GridViewType > BoundaryInfoType;

  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix< RangeFieldType >  MatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector< RangeFieldType >           VectorType;

  DiscretizationBase(const GridPartType& grid_part,
                     const BoundaryInfoType& boundary_info,
                     const FunctionType& diffusion,
                     const FunctionType& force,
                     const FunctionType& dirichlet)
    : grid_part_(grid_part)
    , boundary_info_(boundary_info)
    , diffusion_(diffusion)
    , force_(force)
    , dirichlet_(dirichlet)
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

  template< class ReferenceSolutionType, class DiscreteSolutionType >
  void compute_errors(const ReferenceSolutionType& referencec_solution,
                      const DiscreteSolutionType& discrete_solution) const
  {
    static_assert(std::is_base_of< Dune::Stuff::LocalizableFunction, ReferenceSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    static_assert(std::is_base_of< Dune::Stuff::LocalizableFunction, DiscreteSolutionType >::value,
                  "ReferenceSolutionType has to be derived from Stuff::LocalizableFunction!");
    typedef Dune::GDT::DiscreteFunction::Difference< ReferenceSolutionType, DiscreteSolutionType > DifferenceType;
    DifferenceType difference(referencec_solution, discrete_solution);

  }

private:
  const GridPartType& grid_part_;
  const BoundaryInfoType& boundary_info_;
  const FunctionType& diffusion_;
  const FunctionType& force_;
  const FunctionType& dirichlet_;
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

  CGDiscretization(const GridPartType& grid_part,
                   const BoundaryInfoType& boundary_info,
                   const FunctionType& diffusion,
                   const FunctionType& force,
                   const FunctionType& dirichlet)
    : BaseType(grid_part, boundary_info, diffusion, force, dirichlet)
    , space_(BaseType::grid_part())
  {}

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
  //  //   * L2 neumann functional
  //  typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > L2FaceFunctionalType;
  //  const L2FaceFunctionalType neumannFunctional(neumann);

    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType force_vector(space_.mapper().size());
    auto dirichlet_vector = std::make_shared< VectorType >(space_.mapper().size());
  //  VectorType neumann_nector(space_.mapper().size());
    VectorType rhs_vector(space_.mapper().size());
    auto solution = std::make_shared< DiscreteFunctionType >(space_,
                                                             std::make_shared< VectorType >(space_.mapper().size()));

    // * dirichlet boundary values
    DiscreteFunctionType dirichlet_projection(space_, dirichlet_vector, "dirichlet");

    typedef Operator::DirichletProjection< FunctionType, DiscreteFunctionType > DirichletProjectionOperatorType;
    const DirichletProjectionOperatorType dirichlet_projection_operator(BaseType::boundary_info());
    dirichlet_projection_operator.apply(BaseType::dirichlet(), dirichlet_projection);

    // * local matrix assembler
    typedef LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalEllipticOperatorMatrixAssemblerType;
    const LocalEllipticOperatorMatrixAssemblerType diffusion_matrix_assembler(diffusion_operator);
    // * local vector assemblers
    //   * force vector
    typedef LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalL2VolumeFunctionalVectorAssemblerType;
    const LocalL2VolumeFunctionalVectorAssemblerType force_vector_assembler(force_functional);
  //  //   * neumann vector
  //  typedef LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalL2FaceFunctionalVectorAssemblerType;
  //  const LocalL2FaceFunctionalVectorAssemblerType neumannVectorAssembler(neumannFunctional);
    // * system assembler
    typedef SystemAssembler< SpaceType, SpaceType > SystemAssemblerType;
    SystemAssemblerType system_assembler(space_);
    system_assembler.addLocalAssembler(diffusion_matrix_assembler, system_matrix);
    system_assembler.addLocalAssembler(force_vector_assembler, force_vector);
  //  system_assembler.addLocalAssembler(neumannVectorAssembler,
  //                                    SystemAssemblerType::AssembleOnNeumann(*boundaryInfo),
  //                                    *neumannVector);
    system_assembler.assemble();

    Constraints::Dirichlet< typename GridPartType::GridViewType,
        RangeFieldType > dirichlet_constraints(BaseType::boundary_info(),
                                                                   space_.mapper().maxNumDofs(),
                                                                   space_.mapper().maxNumDofs());
    rhs_vector.backend() = force_vector.backend()
  //      + neumannVector->backend()
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

  void visualize(const VectorType& vector, const std::string filename)
  {
    space_.visualize(vector, filename);
  }

private:
  const SpaceType space_;
};


} // namespace Example

#endif // DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATIONS_HH
