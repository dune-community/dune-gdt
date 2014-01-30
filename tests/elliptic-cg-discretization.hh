// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/common/timer.hh>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/operator/prolongations.hh>

#include "elliptic-testcases.hh"


namespace EllipticCG {


template< class GridPartType, int polynomialOrder >
class Discretization
{
public:
  static const unsigned int             dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype  DomainFieldType;

  static const unsigned int dimRange = 1;
  typedef double            RangeFieldType;

  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::GridboundaryInterface< typename GridPartType::IntersectionType > BoundaryInfoType;
  typedef Dune::Stuff::LocalizableFunctionInterface
      < typename GridPartType::template Codim< 0 >::EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    FunctionType;

  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix< RangeFieldType >  MatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector< RangeFieldType >           VectorType;

  typedef Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper
      < GridPartType, polOrder, RangeFieldType, dimRange > SpaceType;

  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType >      DiscreteFunctionType;
  typedef Dune::GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;

  Discretization(const std::shared_ptr< const GridPartType >& gp,
                 const BoundaryInfoType& info,
                 const FunctionType& diff,
                 const FunctionType& forc,
                 const FunctionType& dir,
                 const FunctionType& neu)
    : space_(gp)
    , boundary_info_(info)
    , diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neu)
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  VectorType create_vector() const
  {
    return VectorType(space_.mapper().size());
  }

  void solve(VectorType& solution) const
  {
    using namespace Dune;
    using namespace Dune::GDT;

    // left hand side
    // * elliptic diffusion operator
    typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< FunctionType > >  EllipticOperatorType;
    const EllipticOperatorType diffusion_operator(diffusion_);
    // right hand side
    // * L2 force functional
    typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > L2VolumeFunctionalType;
    const L2VolumeFunctionalType force_functional(force_);
    // * L2 neumann functional
    typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > L2FaceFunctionalType;
    const L2FaceFunctionalType neumann_functional(neumann_);

    const std::unique_ptr< Stuff::LA::SparsityPatternDefault > sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType dirichlet_vector(space_.mapper().size());
    VectorType rhs_vector(space_.mapper().size());

    // dirichlet boundary values
    DiscreteFunctionType dirichlet_projection(space_, dirichlet_vector, "dirichlet");
    typedef ProjectionOperator::Dirichlet< GridPartType > DirichletProjectionOperatorType;
    const DirichletProjectionOperatorType dirichlet_projection_operator(*(space_.gridPart()), boundary_info_);
    dirichlet_projection_operator.apply(dirichlet_, dirichlet_projection);

    // local matrix assembler
    typedef LocalAssembler::Codim0Matrix< EllipticOperatorType > LocalEllipticOperatorMatrixAssemblerType;
    const LocalEllipticOperatorMatrixAssemblerType diffusion_matrix_assembler(diffusion_operator);
    // local vector assemblers
    // * force vector
    typedef LocalAssembler::Codim0Vector< L2VolumeFunctionalType > LocalL2VolumeFunctionalVectorAssemblerType;
    const LocalL2VolumeFunctionalVectorAssemblerType force_vector_assembler(force_functional);
    // * neumann vector
    typedef LocalAssembler::Codim1Vector< L2FaceFunctionalType > LocalL2FaceFunctionalVectorAssemblerType;
    const LocalL2FaceFunctionalVectorAssemblerType neumann_vector_assembler(neumann_functional);
    // system assembler
    typedef SystemAssembler< SpaceType > SystemAssemblerType;
    SystemAssemblerType system_assembler(space_);
    system_assembler.addLocalAssembler(diffusion_matrix_assembler, system_matrix);
    system_assembler.addLocalAssembler(force_vector_assembler, rhs_vector);
    system_assembler.addLocalAssembler(neumann_vector_assembler,
                                      typename SystemAssemblerType::AssembleOnNeumann(boundary_info_),
                                      rhs_vector);
    system_assembler.assemble();

    Constraints::Dirichlet < typename GridPartType::IntersectionType, RangeFieldType >
      dirichlet_constraints(boundary_info_, space_.mapper().maxNumDofs(), space_.mapper().maxNumDofs());
    rhs_vector.backend() -= system_matrix.backend()*dirichlet_vector.backend();
    system_assembler.addLocalConstraints(dirichlet_constraints, system_matrix);
    system_assembler.addLocalConstraints(dirichlet_constraints, rhs_vector);
    system_assembler.applyConstraints();

    // solve
    Dune::Stuff::LA::Solver< MatrixType > linear_solver(system_matrix);
    const size_t failure = linear_solver.apply(rhs_vector, solution);
    if (failure)
      DUNE_THROW_COLORFULLY(Dune::MathError,
                            "linear solver failed with error code " << failure << " (see dune/stuff/solver.hh)!");
    solution.backend() += dirichlet_vector.backend();
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename, const std::string name) const
  {
    ConstDiscreteFunctionType function(space_, vector, name);
    function.visualize(filename);
  }

private:
  const SpaceType space_;
  const BoundaryInfoType& boundary_info_;
  const FunctionType& diffusion_;
  const FunctionType& force_;
  const FunctionType& dirichlet_;
  const FunctionType& neumann_;
}; // class Discretization


template< class TestCase, int polOrder >
class EocStudy
  : public Dune::Stuff::Common::ConvergenceStudy
{
  typedef Dune::Stuff::Common::ConvergenceStudy BaseType;
protected:
  typedef typename TestCase::GridPartType GridPartType;
  typedef typename TestCase::EntityType   EntityType;

  typedef typename TestCase::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = TestCase::dimDomain;
  typedef typename TestCase::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRange = TestCase::dimRange;

  typedef Discretization< GridPartType, polOrder >  DiscretizationType;
  typedef typename DiscretizationType::VectorType                 VectorType;
  typedef typename DiscretizationType::DiscreteFunctionType       DiscreteFunctionType;
  typedef typename DiscretizationType::ConstDiscreteFunctionType  ConstDiscreteFunctionType;

  typedef typename TestCase::ExactSolutionType ExactSolutionType;

public:
  EocStudy(const TestCase& test)
    : test_(test)
    , current_level_(0)
    , last_computed_level_(666)
    , reference_solution_computed_(false)
    , reference_discretization_(nullptr)
    , reference_solution_vector_(nullptr)
    , current_solution_vector_(nullptr)
  {}

  virtual ~EocStudy() {}

  virtual std::string identifier() const DS_OVERRIDE
  {
    return "continuous galerkin discretization, polOrder " + Dune::Stuff::Common::toString(polOrder);
  }

  virtual size_t num_refinements() const DS_OVERRIDE
  {
    if (test_.num_levels() == 0)
      return test_.num_levels();
    else
      return test_.num_levels() - 1;
  }

  virtual std::vector< std::string > provided_norms() const DS_OVERRIDE
  {
    return {"L2", "H1_semi"};
  }

  virtual size_t expected_rate(const std::string type) const
  {
    if (type.compare("L2") == 0)
      return polOrder + 1;
    else if (type.compare("H1_semi") == 0)
      return polOrder;
    else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  virtual double norm_reference_solution(const std::string type) DS_OVERRIDE
  {
    if (test_.provides_exact_solution()) {
      return compute_norm(*(test_.reference_grid_part()), test_.exact_solution(), type);
    } else {
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(reference_discretization_->space(),
                                                         *reference_solution_vector_,
                                                         "CG reference solution");
      // compute norm
      return compute_norm(*(test_.reference_grid_part()), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const DS_OVERRIDE
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return grid_part->grid().size(grid_part->level(), 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const DS_OVERRIDE
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return Dune::Fem::GridWidth::calcGridWidth(*grid_part);
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() DS_OVERRIDE
  {
    using namespace Dune;
    using namespace Dune::GDT;
    Timer timer;
    double elapsed = 0;
    if (current_level_ != last_computed_level_) {
      assert(current_level_ < test_.num_levels());
      // compute solution
      const auto grid_part = test_.level_grid_part(current_level_);
      const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                              test_.dirichlet(), test_.neumann());
      auto current_level_solution_vector = discretization.create_vector();
      discretization.solve(current_level_solution_vector);
      const ConstDiscreteFunctionType current_level_solution(discretization.space(),
                                                             current_level_solution_vector,
                                                             "solution on current level");
      // prolong to reference grid part
      elapsed += timer.elapsed();
      if (!reference_solution_computed_)
        compute_reference_solution();
      timer.reset();
      const auto reference_grid_part = test_.reference_grid_part();
      const ProlongationOperator::Generic< GridPartType > prolongation_operator(*reference_grid_part);
      assert(reference_discretization_);
      current_solution_vector_
          = std::unique_ptr< VectorType >(new VectorType(reference_discretization_->create_vector()));
      DiscreteFunctionType reference_level_solution(reference_discretization_->space(),
                                                    *current_solution_vector_,
                                                    "solution on reference grid part");
      prolongation_operator.apply(current_level_solution, reference_level_solution);
      last_computed_level_ = current_level_;
    }
    return timer.elapsed() + elapsed;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) DS_OVERRIDE
  {
    // get current solution
    assert(current_level_ < test_.num_levels());
    if (last_computed_level_ != current_level_) {
      compute_on_current_refinement();
    }
    assert(last_computed_level_ == current_level_);
    assert(current_solution_vector_);
    if (!reference_solution_computed_)
      compute_reference_solution();
    assert(reference_discretization_);
    const ConstDiscreteFunctionType current_solution(reference_discretization_->space(),
                                                     *current_solution_vector_,
                                                     "current solution");
    // compute error
    if (test_.provides_exact_solution()) {
      typedef Dune::Stuff::Function::Difference< ExactSolutionType, ConstDiscreteFunctionType > DifferenceType;
      const DifferenceType difference(test_.exact_solution(), current_solution);
      return compute_norm(*(test_.reference_grid_part()), difference, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(reference_discretization_->space(),
                                                         *reference_solution_vector_,
                                                         "CG reference solution");
      typedef Dune::Stuff::Function::Difference< ConstDiscreteFunctionType, ConstDiscreteFunctionType > DifferenceType;
      const DifferenceType difference(reference_solution, current_solution);
      return compute_norm(*(test_.reference_grid_part()), difference, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() DS_OVERRIDE
  {
    if (current_level_ < test_.num_levels())
      ++current_level_;
  }

  std::vector< double > expected_results(const std::string type) const
  {
    if (std::is_same< TestCase, EllipticTestCase::ESV07< Dune::ALUConformGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.97e-01, 4.86e-02, 1.22e-02, 3.03e-03};
        else if (type.compare("H1_semi") == 0)
          return {4.12e-01, 2.08e-01, 1.04e-01, 5.18e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::LocalThermalBlock< Dune::ALUConformGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {7.93e-02, 4.16e-02, 1.20e-02, 2.73e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.52e-01, 3.03e-01, 1.64e-01, 7.52e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::ER07< Dune::ALUConformGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.49e-01, 3.83e-02, 9.66e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.60e-01, 1.85e-01, 9.25e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::MixedBoundaryTypes< Dune::ALUConformGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {8.32e-02, 2.23e-02, 5.53e-03, 1.20e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.12e-01, 1.65e-01, 8.24e-02, 3.76e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {2.91e-03, 1.13e-03, 3.72e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.37e-01, 1.43e-01, 7.57e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  } // ... expected_results(...)

  virtual std::map< std::string, std::vector< double > > run(std::ostream& out = std::cout)
  {
    return BaseType::run(true, out);
  }

private:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_
          = std::unique_ptr<DiscretizationType>(new DiscretizationType(test_.reference_grid_part(),
                                                                       test_.boundary_info(),
                                                                       test_.diffusion(),
                                                                       test_.force(),
                                                                       test_.dirichlet(),
                                                                       test_.neumann()));
      reference_solution_vector_
          = std::unique_ptr< VectorType >(new VectorType(reference_discretization_->create_vector()));
      reference_discretization_->solve(*reference_solution_vector_);
      reference_solution_computed_ = true;
    }
  } // ... compute_reference_solution()

  template< class GridPartType, class FunctionType >
  double compute_norm(const GridPartType& grid_part, const FunctionType& function, const std::string type)
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (type.compare("L2") == 0) {
      ProductOperator::L2< GridPartType > l2_product_operator(grid_part);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      ProductOperator::H1Semi< GridPartType > h1_product_operator(grid_part);
      return std::sqrt(h1_product_operator.apply2(function, function));
    } else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

  const TestCase& test_;
  size_t current_level_;
  size_t last_computed_level_;
  bool reference_solution_computed_;
  std::unique_ptr< DiscretizationType > reference_discretization_;
  std::unique_ptr< VectorType > reference_solution_vector_;
  std::unique_ptr< VectorType > current_solution_vector_;
}; // class EOCStudy


} // namespace EllipticCG

#endif // DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH
