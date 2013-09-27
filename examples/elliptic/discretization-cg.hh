// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_CG_HH
#define DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_CG_HH

#include "discretization-base.hh"

namespace Discretization {


template< class GridPartType, int polOrder >
class ContinuousGalerkin
  : public Base< GridPartType, polOrder >
{
  typedef Base< GridPartType, polOrder > BaseType;
public:
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = BaseType::dimRange;
  typedef typename BaseType::FunctionType   FunctionType;

  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::VectorType VectorType;

  typedef Dune::GDT::ContinuousLagrangeSpace::FemWrapper< GridPartType, polOrder, RangeFieldType, dimRange >  SpaceType;
  typedef Dune::GDT::DiscreteFunctionDefault< SpaceType, VectorType > DiscreteFunctionType;

  ContinuousGalerkin(const GridPartType& gp,
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


} // namespace Discretization

#endif // DUNE_GDT_EXAMPLES_ELLIPTIC_DISCRETIZATION_CG_HH
