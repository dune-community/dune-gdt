// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_SWIPDG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_SWIPDG_DISCRETIZATION_HH

#include <memory>
#include <vector>
#include <type_traits>
#include <math.h>
#include <limits>

#include <dune/common/timer.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver/eigen.hh>
#include <dune/stuff/la/solver/fasp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/raviartthomas/fem-localfunctions.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localevaluation/swipdg.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/operator/products.hh>
#include <dune/gdt/operator/prolongations.hh>
#include <dune/gdt/operator/reconstructions.hh>

#include "elliptic-testcases.hh"


namespace EllipticSWIPDG {


template <class GridPartType, int polynomialOrder>
class Discretization
{
public:
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype DomainFieldType;

  static const unsigned int dimRange = 1;
  typedef double RangeFieldType;

  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::GridboundaryInterface<typename GridPartType::IntersectionType> BoundaryInfoType;
  typedef Dune::Stuff::LocalizableFunctionInterface<typename GridPartType::template Codim<0>::EntityType,
                                                    DomainFieldType, dimDomain, RangeFieldType, dimRange> FunctionType;

  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix<RangeFieldType> MatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector<RangeFieldType> VectorType;

  typedef Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, polOrder, RangeFieldType,
                                                                          dimRange> SpaceType;

  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;

  Discretization(const std::shared_ptr<const GridPartType>& gp, const BoundaryInfoType& info, const FunctionType& diff,
                 const FunctionType& forc, const FunctionType& dir, const FunctionType& neu)
    : space_(gp)
    , boundary_info_(info)
    , diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neu)
    , beta_(1.0)
  {
  }

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

    const std::unique_ptr<Stuff::LA::SparsityPatternDefault> sparsity_pattern(space_.computePattern());
    MatrixType system_matrix(space_.mapper().size(), space_.mapper().size(), *sparsity_pattern);
    VectorType rhs_vector(space_.mapper().size());
    typedef SystemAssembler<SpaceType> SystemAssemblerType;
    SystemAssemblerType systemAssembler(space_);

    // volume terms
    // * lhs
    typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<FunctionType>> EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(diffusion_);
    const LocalAssembler::Codim0Matrix<EllipticOperatorType> diffusionMatrixAssembler(ellipticOperator);
    systemAssembler.addLocalAssembler(diffusionMatrixAssembler, system_matrix);
    // * rhs
    typedef LocalFunctional::Codim0Integral<LocalEvaluation::Product<FunctionType>> ForceFunctionalType;
    const ForceFunctionalType forceFunctional(force_);
    const LocalAssembler::Codim0Vector<ForceFunctionalType> forceVectorAssembler(forceFunctional);
    systemAssembler.addLocalAssembler(forceVectorAssembler, rhs_vector);
    // inner face terms
    typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDG::Inner<FunctionType>> CouplingOperatorType;
    const CouplingOperatorType couplingOperator(diffusion_, beta_);
    const LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> couplingMatrixAssembler(couplingOperator);
    systemAssembler.addLocalAssembler(
        couplingMatrixAssembler, typename SystemAssemblerType::AssembleOnInnerPrimally(), system_matrix);
    // dirichlet boundary face terms
    // * lhs
    typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SWIPDG::BoundaryLHS<FunctionType>>
        DirichletOperatorType;
    const DirichletOperatorType dirichletOperator(diffusion_, beta_);
    const LocalAssembler::Codim1BoundaryMatrix<DirichletOperatorType> dirichletMatrixAssembler(dirichletOperator);
    systemAssembler.addLocalAssembler(
        dirichletMatrixAssembler, typename SystemAssemblerType::AssembleOnDirichlet(boundary_info_), system_matrix);
    // * rhs
    typedef LocalFunctional::Codim1Integral<LocalEvaluation::SWIPDG::BoundaryRHS<FunctionType, FunctionType>>
        DirichletFunctionalType;
    const DirichletFunctionalType dirichletFunctional(diffusion_, dirichlet_, beta_);
    const LocalAssembler::Codim1Vector<DirichletFunctionalType> dirichletVectorAssembler(dirichletFunctional);
    systemAssembler.addLocalAssembler(
        dirichletVectorAssembler, typename SystemAssemblerType::AssembleOnDirichlet(boundary_info_), rhs_vector);
    // neumann boundary face terms
    // * rhs
    typedef LocalFunctional::Codim1Integral<LocalEvaluation::Product<FunctionType>> NeumannFunctionalType;
    const NeumannFunctionalType neumannFunctional(neumann_);
    const LocalAssembler::Codim1Vector<NeumannFunctionalType> neumannVectorAssembler(neumannFunctional);
    systemAssembler.addLocalAssembler(
        neumannVectorAssembler, typename SystemAssemblerType::AssembleOnNeumann(boundary_info_), rhs_vector);
    // do all the work
    systemAssembler.assemble();

    // solve
    std::unique_ptr<typename Dune::Stuff::LA::SolverInterface<MatrixType, VectorType>> linear_solver(nullptr);
    Dune::ParameterTree linear_solver_settings;
#ifdef HAVE_FASP
    typedef typename Dune::Stuff::LA::AmgSolver<MatrixType, VectorType> LinearSolverType;
    linear_solver_settings              = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-10";
    linear_solver_settings["maxIter"]   = Dune::Stuff::Common::toString(space_.mapper().size());
    linear_solver                       = std::unique_ptr<LinearSolverType>(new LinearSolverType());
#else
    typedef typename Dune::Stuff::LA::BicgstabILUTSolver<MatrixType, VectorType> LinearSolverType;
    linear_solver_settings              = LinearSolverType::defaultSettings();
    linear_solver_settings["precision"] = "1e-10";
    linear_solver_settings["maxIter"]   = Dune::Stuff::Common::toString(space_.mapper().size());
    linear_solver                       = std::unique_ptr<LinearSolverType>(new LinearSolverType());
#endif
    assert(linear_solver);
    const size_t failure = linear_solver->apply(system_matrix, rhs_vector, solution, linear_solver_settings);
    if (failure)
      DUNE_THROW(Dune::MathError, "\nERROR: linear solver reported a problem!");
    if (solution.size() != space_.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver produced a solution of wrong size (is " << solution.size() << ", should be "
                                                                                 << space_.mapper().size()
                                                                                 << ")!");
  } // ... solve()

  void visualize(const VectorType& vector, const std::string filename, const std::string name) const
  {
    ConstDiscreteFunctionType function(space_, vector, name);
    function.visualize(space_.gridPart()->gridView(), filename);
  }

private:
  const SpaceType space_;
  const BoundaryInfoType& boundary_info_;
  const FunctionType& diffusion_;
  const FunctionType& force_;
  const FunctionType& dirichlet_;
  const FunctionType& neumann_;
  const RangeFieldType beta_;
}; // class Discretization


template <class TestCase, int polOrder>
class EocStudy : public Dune::Stuff::Common::ConvergenceStudy
{
protected:
  typedef Dune::Stuff::Common::ConvergenceStudy BaseType;

  typedef typename TestCase::GridPartType GridPartType;
  typedef typename TestCase::EntityType EntityType;

  typedef typename TestCase::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = TestCase::dimDomain;
  typedef typename TestCase::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCase::dimRange;

  typedef Discretization<GridPartType, polOrder> DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscretizationType::ConstDiscreteFunctionType ConstDiscreteFunctionType;

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
  {
  }

  virtual ~EocStudy()
  {
  }

  virtual std::string identifier() const override
  {
    return "SWIP discontinuous galerkin discretization, polOrder " + Dune::Stuff::Common::toString(polOrder);
  }

  virtual size_t num_refinements() const override
  {
    if (test_.num_levels() == 0)
      return test_.num_levels();
    else
      return test_.num_levels() - 1;
  }

  virtual std::vector<std::string> provided_norms() const override
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

  virtual double norm_reference_solution(const std::string type) override
  {
    if (test_.provides_exact_solution()) {
      return compute_norm(*(test_.reference_grid_part()), test_.exact_solution(), type);
    } else {
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "CG reference solution");
      // compute norm
      return compute_norm(*(test_.reference_grid_part()), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const override
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return grid_part->grid().size(grid_part->level(), 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const override
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return Dune::Fem::GridWidth::calcGridWidth(*grid_part);
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() override
  {
    using namespace Dune;
    using namespace Dune::GDT;
    Timer timer;
    double elapsed = 0;
    if (current_level_ != last_computed_level_) {
      assert(current_level_ < test_.num_levels());
      // compute solution
      const auto grid_part = test_.level_grid_part(current_level_);
      const DiscretizationType discretization(
          grid_part, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
      current_solution_vector_on_level_ = std::unique_ptr<VectorType>(new VectorType(discretization.create_vector()));
      discretization.solve(*current_solution_vector_on_level_);
      const ConstDiscreteFunctionType current_level_solution(
          discretization.space(), *current_solution_vector_on_level_, "solution on current level");
      // prolong to reference grid part
      elapsed += timer.elapsed();
      if (!reference_solution_computed_)
        compute_reference_solution();
      timer.reset();
      const auto reference_grid_part = test_.reference_grid_part();
      const ProlongationOperator::L2<GridPartType> prolongation_operator(*reference_grid_part);
      assert(reference_discretization_);
      current_solution_vector_ =
          std::unique_ptr<VectorType>(new VectorType(reference_discretization_->create_vector()));
      DiscreteFunctionType reference_level_solution(
          reference_discretization_->space(), *current_solution_vector_, "solution on reference grid part");
      prolongation_operator.apply(current_level_solution, reference_level_solution);
      last_computed_level_ = current_level_;
    }
    return timer.elapsed() + elapsed;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override
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
    const ConstDiscreteFunctionType current_solution(
        reference_discretization_->space(), *current_solution_vector_, "current solution");
    // compute error
    if (test_.provides_exact_solution()) {
      typedef Dune::Stuff::Function::Difference<ExactSolutionType, ConstDiscreteFunctionType> DifferenceType;
      const DifferenceType difference(test_.exact_solution(), current_solution);
      return compute_norm(*(test_.reference_grid_part()), difference, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "CG reference solution");
      typedef Dune::Stuff::Function::Difference<ConstDiscreteFunctionType, ConstDiscreteFunctionType> DifferenceType;
      const DifferenceType difference(reference_solution, current_solution);
      return compute_norm(*(test_.reference_grid_part()), difference, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override
  {
    if (current_level_ < test_.num_levels())
      ++current_level_;
  }

  std::vector<double> expected_results(const std::string type) const
  {
    if (std::is_same<TestCase, EllipticTestCase::ESV07<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.15e-01, 3.04e-02, 7.51e-03, 1.86e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.79e-01, 1.90e-01, 9.38e-02, 4.67e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {1.25e-02, 1.42e-03, 1.69e-04, 2.08e-05};
        else if (type.compare("H1_semi") == 0)
          return {7.84e-02, 2.01e-02, 5.02e-03, 1.26e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::LocalThermalBlock<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {5.85e-02, 2.12e-02, 5.89e-03, 1.38e-03};
        else if (type.compare("H1_semi") == 0)
          return {4.49e-01, 3.06e-01, 1.58e-01, 6.83e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {1.27e-02, 2.26e-03, 4.15e-04, 8.21e-05};
        else if (type.compare("H1_semi") == 0)
          return {1.78e-01, 6.26e-02, 2.05e-02, 6.35e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::ER07<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {6.10e-02, 1.66e-02, 4.23e-03};
        else if (type.compare("H1_semi") == 0)
          return {2.99e-01, 1.47e-01, 7.24e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {6.43e-03, 8.24e-04, 1.05e-04};
        else if (type.compare("H1_semi") == 0)
          return {5.41e-02, 1.42e-02, 3.56e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::MixedBoundaryTypes<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {5.21e-02, 1.45e-02, 3.65e-03, 8.15e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.76e-01, 1.42e-01, 6.92e-02, 3.09e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {4.64e-03, 7.69e-04, 1.38e-04, 2.84e-05};
        else if (type.compare("H1_semi") == 0)
          return {4.56e-02, 1.57e-02, 5.78e-03, 2.20e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::Spe10Model1<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.17e-03, 3.69e-04};
        else if (type.compare("H1_semi") == 0)
          return {1.66e-01, 8.26e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {4.05e-04, 1.46e-04};
        else if (type.compare("H1_semi") == 0)
          return {7.80e-02, 3.67e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  }

  virtual std::map<std::string, std::vector<double>> run(std::ostream& out = std::cout)
  {
    return BaseType::run(true, out);
  }

private:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_ =
          std::unique_ptr<DiscretizationType>(new DiscretizationType(test_.reference_grid_part(),
                                                                     test_.boundary_info(),
                                                                     test_.diffusion(),
                                                                     test_.force(),
                                                                     test_.dirichlet(),
                                                                     test_.neumann()));
      reference_solution_vector_ =
          std::unique_ptr<VectorType>(new VectorType(reference_discretization_->create_vector()));
      reference_discretization_->solve(*reference_solution_vector_);
      reference_solution_computed_ = true;
    }
  } // ... compute_reference_solution()

  template <class GridPartType, class FunctionType>
  double compute_norm(const GridPartType& grid_part, const FunctionType& function, const std::string type)
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (type.compare("L2") == 0) {
      ProductOperator::L2<GridPartType> l2_product_operator(grid_part);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      ProductOperator::H1Semi<GridPartType> h1_product_operator(grid_part);
      return std::sqrt(h1_product_operator.apply2(function, function));
    } else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

protected:
  const TestCase& test_;
  size_t current_level_;
  size_t last_computed_level_;
  bool reference_solution_computed_;
  std::unique_ptr<DiscretizationType> reference_discretization_;
  std::unique_ptr<VectorType> reference_solution_vector_;
  std::unique_ptr<VectorType> current_solution_vector_on_level_;
  std::unique_ptr<VectorType> current_solution_vector_;
}; // class EOCStudy


template <class TestCase>
class EstimatorStudy : public EocStudy<TestCase, 1>
{
  typedef EocStudy<TestCase, 1> BaseType;
  static const unsigned int polOrder = 1;

  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = BaseType::dimRange;

  typedef typename BaseType::DiscretizationType DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscretizationType::ConstDiscreteFunctionType ConstDiscreteFunctionType;

  typedef typename TestCase::ExactSolutionType ExactSolutionType;

  static std::string residual_estimator_ESV07_string()
  {
    return "eta_R (ESV07)";
  }
  static std::string residual_estimator_ESV10_string()
  {
    return "eta_R (ESV10)";
  }
  static std::string diffusive_flux_estimator_ESV10_string()
  {
    return "eta_DF (ESV10)";
  }
  static std::string nonconformity_estimator_ESV10_string()
  {
    return "eta_NC (ESV10)";
  }
  static std::string estimator_ESV10_string()
  {
    return "estimator (ESV10)";
  }

public:
  EstimatorStudy(const TestCase& test)
    : BaseType(test)
  {
  }

  virtual ~EstimatorStudy()
  {
  }

  virtual std::vector<std::string> provided_norms() const override
  {
    return {/*"L2",*/ "H1_semi",
            /*residual_estimator_ESV07_string(),*/ residual_estimator_ESV10_string(),
            diffusive_flux_estimator_ESV10_string(),
            nonconformity_estimator_ESV10_string(),
            estimator_ESV10_string()};
  }

  virtual size_t expected_rate(const std::string type) const
  {
    if (type.compare(residual_estimator_ESV07_string()) == 0)
      return polOrder + 1;
    else if (type.compare(residual_estimator_ESV10_string()) == 0)
      return polOrder + 1;
    else if (type.compare(diffusive_flux_estimator_ESV10_string()) == 0)
      return polOrder;
    else if (type.compare(nonconformity_estimator_ESV10_string()) == 0)
      return polOrder;
    else if (type.compare(estimator_ESV10_string()) == 0)
      return polOrder + 1;
    else
      return BaseType::expected_rate(type);
  } // ... expected_rate(...)

  virtual double current_error_norm(const std::string type) override
  {
    if (type.compare(residual_estimator_ESV07_string()) == 0) {
      return compute_residual_estimator_ESV07();
    } else if (type.compare(residual_estimator_ESV10_string()) == 0) {
      return compute_residual_estimator_ESV10();
    } else if (type.compare(diffusive_flux_estimator_ESV10_string()) == 0) {
      return compute_diffusive_flux_estimator_ESV10();
    } else if (type.compare(nonconformity_estimator_ESV10_string()) == 0) {
      return compute_nonconformity_estimator_ESV10();
    } else if (type.compare(estimator_ESV10_string()) == 0) {
      return compute_estimator_ESV10();
    } else
      return BaseType::current_error_norm(type);
  } // ... current_error_norm(...)

  std::vector<double> expected_results(const std::string type) const
  {
    return {0.0};
    //    if (std::is_same< TestCase, EllipticTestCase::ESV07< Dune::ALUConformGrid< 2, 2 > > >::value) {
    //      if (polOrder == 1) {
    //        if (type.compare("L2") == 0)
    //          return {1.97e-01, 4.86e-02, 1.22e-02, 3.03e-03};
    //        else if (type.compare("H1") == 0)
    //          return {4.12e-01, 2.08e-01, 1.04e-01, 5.18e-02};
    //        else
    //          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
    //      } else
    //        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    //    } else if (std::is_same< TestCase, EllipticTestCase::LocalThermalBlock< Dune::ALUConformGrid< 2, 2 > >
    //    >::value) {
    //      if (polOrder == 1) {
    //        if (type.compare("L2") == 0)
    //          return {8.04e-02, 4.39e-02, 1.27e-02, 2.88e-03};
    //        else if (type.compare("H1") == 0)
    //          return {3.61e-01, 3.15e-01, 1.71e-01, 7.83e-02};
    //        else
    //          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
    //      } else
    //        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    //    } else if (std::is_same< TestCase, EllipticTestCase::ER07< Dune::ALUConformGrid< 2, 2 > > >::value) {
    //      if (polOrder == 1) {
    //        if (type.compare("L2") == 0)
    //          return {1.93e-01, 5.73e-02, 1.55e-02};
    //        else if (type.compare("H1") == 0)
    //          return {3.62e-01, 1.85e-01, 9.25e-02};
    //        else
    //          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
    //      } else
    //        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    //    } else if (std::is_same< TestCase, EllipticTestCase::MixedBoundaryTypes< Dune::ALUConformGrid< 2, 2 > >
    //    >::value) {
    //      if (polOrder == 1) {
    //        if (type.compare("L2") == 0)
    //          return {1.19e-01, 3.12e-02, 7.73e-03, 1.66e-03};
    //        else if (type.compare("H1") == 0)
    //          return {3.21e-01, 1.67e-01, 8.30e-02, 3.76e-02};
    //        else
    //          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
    //      } else
    //        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    //    } else if (std::is_same< TestCase, EllipticTestCase::Spe10Model1< Dune::ALUConformGrid< 2, 2 > > >::value) {
    //      if (polOrder == 1) {
    //        if (type.compare("L2") == 0)
    //          return {2.91e-03, 1.13e-03, 3.72e-04};
    //        else if (type.compare("H1") == 0)
    //          return {2.37e-01, 1.43e-01, 7.57e-02};
    //        else
    //          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
    //      } else
    //        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    //    } else
    //      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType
    //      combination!");
  }

  virtual std::map<std::string, std::vector<double>> run(std::ostream& out = std::cout)
  {
    return BaseType::BaseType::run(false, out);
  }

private:
  double compute_residual_estimator_ESV07()
  {
    using namespace Dune;
    typename TestCase::DiffusionType::RangeType diffusion_min = std::numeric_limits<RangeFieldType>::max();
    long double ret                                           = 0;
    // walk the grid
    const auto grid_part     = test_.level_grid_part(current_level_);
    const auto entity_it_end = grid_part->template end<0>();
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get local functions
      const auto local_force      = test_.force().local_function(entity);
      const auto local_diffusion  = test_.diffusion().local_function(entity);
      const auto force_mean_value = local_force->evaluate(entity.geometry().local(entity.geometry().center()));
      diffusion_min               = std::numeric_limits<RangeFieldType>::max();
      // do a volume quadrature
      const size_t quadrature_order = 2 * local_force->order() + 4;
      const auto& quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * quadrature_order + 1);
      double integral = 0;
      for (auto quadrature_point : quadrature) {
        const auto& point_entity        = quadrature_point.position();
        const double quadrature_weight  = quadrature_point.weight();
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        // evaluate
        const auto diffusion_value = local_diffusion->evaluate(point_entity);
        if (diffusion_value < diffusion_min)
          diffusion_min        = diffusion_value;
        const auto force_value = local_force->evaluate(point_entity);
        const auto difference = force_value - force_mean_value;
        integral += quadrature_weight * integration_factor * (difference * difference);
      } // do a volume quadrature
      // apply scaling factor
      ret += (std::pow(entity.geometry().volume(), 2) / (diffusion_min * M_PIl * M_PIl)) * integral;
    } // walk the grid
    return std::sqrt(ret);
  } // ... compute_residual_estimator_ESV07(...)

  double compute_residual_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const RangeFieldType poincare_constant_sqrt = RangeFieldType(1) / M_PIl;
    const size_t integration_order              = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(
        grid_part, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(
        discretization.space(), *current_solution_vector_on_level_, "discrete solution");

    typedef RaviartThomasSpace::FemLocalfunctionsWrapper<GridPartType, 0, RangeFieldType, dimDomain> RTN_SpaceType;
    const RTN_SpaceType rtn_space(grid_part);
    typedef DiscreteFunction<RTN_SpaceType, VectorType> RTN_DiscreteFunctionType;
    VectorType rtn_vector(rtn_space.mapper().size());
    RTN_DiscreteFunctionType diffusive_flux_reconstruction(rtn_space, rtn_vector, "diffusive flux reconstruction");
    std::vector<typename RTN_SpaceType::BaseFunctionSetType::RangeType> rtn_basis_values(
        rtn_space.mapper().maxNumDofs(), typename RTN_SpaceType::BaseFunctionSetType::RangeType(0));

    // walk the grid for the second time
    const auto entity_it_end = grid_part->template end<0>();
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_solution_entity     = discrete_solution.local_function(entity);
      const auto local_diffusion_entity    = test_.diffusion().local_function(entity);
      auto local_diffusive_flux            = diffusive_flux_reconstruction.local_discrete_function(entity);
      auto local_diffusive_flux_DoF_vector = local_diffusive_flux.vector();
      const auto rtn_basis                 = rtn_space.baseFunctionSet(entity);
      // get the local finite elements
      const auto rtn_finite_element      = rtn_space.backend().finiteElement(entity);
      const auto& rtn_local_coefficients = rtn_finite_element.localCoefficients();
      std::vector<size_t> intersection_id_to_local_DoF_id_map(rtn_local_coefficients.size(), 0);
      for (size_t ii = 0; ii < rtn_local_coefficients.size(); ++ii)
        intersection_id_to_local_DoF_id_map[rtn_local_coefficients.localKey(ii).subEntity()] = ii;
      // to compute the diffusive flux reconstruction
      // * loop over all intersections
      const auto intersection_end_it = grid_part->iend(entity);
      for (auto intersection_it = grid_part->ibegin(entity); intersection_it != intersection_end_it;
           ++intersection_it) {
        // get intersection information
        const auto& intersection        = *intersection_it;
        const size_t intersection_index = intersection.indexInInside();
        assert(intersection_index < intersection_id_to_local_DoF_id_map.size() && "This should not happen!");
        const size_t intersection_DoF_index = intersection_id_to_local_DoF_id_map[intersection_index];
        // prepare quadrature
        const size_t quadrature_order = 1;
        const auto& face_quadrature =
            QuadratureRules<DomainFieldType, dimDomain - 1>::rule(intersection.type(), 2 * quadrature_order + 1);
        RangeFieldType lhs_integral = 0;
        RangeFieldType rhs_integral = 0;

        if (intersection.boundary() && !intersection.neighbor()) {
          // we are on the domain boundary
          // loop over all quadrature points
          for (auto quadrature_point : face_quadrature) {
            const FieldVector<DomainFieldType, dimDomain - 1>& point_intersection = quadrature_point.position();
            const FieldVector<DomainFieldType, dimDomain> point_entity =
                intersection.geometryInInside().global(point_intersection);
            const RangeFieldType integration_factor = intersection.geometry().integrationElement(point_intersection);
            const RangeFieldType quadrature_weight  = quadrature_point.weight();
            // evaluate
            const auto unit_outer_normal = intersection.unitOuterNormal(point_intersection);
            rtn_basis.evaluate(point_entity, rtn_basis_values);
            const auto& rtn_basis_value  = rtn_basis_values[intersection_DoF_index];
            const auto solution_value    = local_solution_entity->evaluate(point_entity);
            const auto solution_gradient = local_solution_entity->jacobian(point_entity);
            const auto diffusion_value   = local_diffusion_entity->evaluate(point_entity);
            const RangeFieldType sigma   = 14.0;
            const RangeFieldType gamma   = diffusion_value;
            const RangeFieldType penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), 1.0);
            // compute integrals
            lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
            rhs_integral +=
                integration_factor * quadrature_weight
                * (-1.0 * diffusion_value * (solution_gradient[0] * unit_outer_normal) + penalty * solution_value);
          } // loop over all quadrature points
          // set local DoF
          local_diffusive_flux_DoF_vector.set(intersection_DoF_index, rhs_integral / lhs_integral);
        } else if (intersection.neighbor() && !intersection.boundary()) {
          // we are on the inside
          // get the neighbour
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour    = *neighbour_ptr;
          // work only once on each intersection
          if (grid_part->indexSet().index(entity) < grid_part->indexSet().index(neighbour)) {
            // get the neighbours local functions
            const auto local_solution_neighbour  = discrete_solution.local_function(neighbour);
            const auto local_diffusion_neighbour = test_.diffusion().local_function(neighbour);
            // loop over all quadrature points
            for (auto quadrature_point : face_quadrature) {
              const FieldVector<DomainFieldType, dimDomain - 1> point_intersection = quadrature_point.position();
              const FieldVector<DomainFieldType, dimDomain> point_entity =
                  intersection.geometryInInside().global(point_intersection);
              const FieldVector<DomainFieldType, dimDomain> point_neighbour =
                  intersection.geometryInOutside().global(point_intersection);
              const RangeFieldType integration_factor = intersection.geometry().integrationElement(point_intersection);
              const RangeFieldType quadrature_weight  = quadrature_point.weight();
              // evaluate
              const auto unit_outer_normal = intersection.unitOuterNormal(point_intersection);
              rtn_basis.evaluate(point_entity, rtn_basis_values);
              const auto& rtn_basis_value            = rtn_basis_values[intersection_DoF_index];
              const auto solution_value_entity       = local_solution_entity->evaluate(point_entity);
              const auto solution_value_neighbour    = local_solution_entity->evaluate(point_neighbour);
              const auto solution_gradient_entity    = local_solution_entity->jacobian(point_entity);
              const auto solution_gradient_neighbour = local_solution_neighbour->jacobian(point_neighbour);
              const auto diffusion_value_entity      = local_diffusion_entity->evaluate(point_entity);
              const auto diffusion_value_neighbour   = local_diffusion_neighbour->evaluate(point_neighbour);
              // evaluate penalty parameter and weighting
              const RangeFieldType sigma        = 8.0;
              const RangeFieldType delta_plus   = diffusion_value_neighbour;
              const RangeFieldType delta_minus  = diffusion_value_entity;
              const RangeFieldType gamma        = (delta_plus * delta_minus) / (delta_plus + delta_minus);
              const RangeFieldType penalty      = (sigma * gamma) / std::pow(intersection.geometry().volume(), 1.0);
              const RangeFieldType weight_plus  = delta_minus / (delta_plus + delta_minus);
              const RangeFieldType weight_minus = delta_plus / (delta_plus + delta_minus);
              // compute integrals
              lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
              rhs_integral +=
                  integration_factor * quadrature_weight
                  * (-weight_minus * (diffusion_value_entity * (solution_gradient_entity[0] * unit_outer_normal))
                     - weight_plus * (diffusion_value_neighbour * (solution_gradient_neighbour[0] * unit_outer_normal))
                     + penalty * (solution_value_entity - solution_value_neighbour));
            } // loop over all quadrature points
            // set local DoF
            local_diffusive_flux_DoF_vector.set(intersection_DoF_index, rhs_integral / lhs_integral);
          } // work only once on each intersection
        } else
          DUNE_THROW(InvalidStateException, "unknwon intersection type encountered!");
      } // loop over all intersections
    } // walk the grid for the second time

    VectorType diffusive_flux_reconstruction_vector_2(rtn_space.mapper().size());
    DiscreteFunction<RTN_SpaceType, VectorType> diffusive_flux_reconstruction_2(
        rtn_space, diffusive_flux_reconstruction_vector_2, "diffusive flux reconstruction");
    // reconstruct
    typedef ReconstructionOperator::DiffusiveFlux<GridPartType, typename TestCase::DiffusionType>
        ReconstructionOperatorType;
    const ReconstructionOperatorType reconstruction_operator(*grid_part, test_.diffusion());
    reconstruction_operator.apply(discrete_solution, diffusive_flux_reconstruction_2);
    //    const auto difference = diffusive_flux_reconstruction.vector().backend() -
    //    diffusive_flux_reconstruction_2.vector().backend();
    //    const double max = std::max(std::abs(difference.maxCoeff()), std::abs(difference.minCoeff()));
    //    if (max > 0.0)
    //      DUNE_THROW(InvalidStateException, "boom " << max);

    // walk the grid for the third time
    std::vector<RangeFieldType> estimators_residual(grid_part->indexSet().size(0), RangeFieldType(0));
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion                     = test_.diffusion().local_function(entity);
      const auto local_force                         = test_.force().local_function(entity);
      const auto local_diffusive_flux_reconstruction = diffusive_flux_reconstruction_2.local_function(entity);
      RangeFieldType diffusion_min                   = std::numeric_limits<RangeFieldType>::max();
      // do a volume quadrature
      const auto& volume_quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector<DomainFieldType, dimDomain>& point_entity = quadrature_point.position();
        const double quadrature_weight                              = quadrature_point.weight();
        // evaluate
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value = local_diffusion->evaluate(point_entity);
        if (diffusion_value < diffusion_min)
          diffusion_min                    = diffusion_value;
        const auto force_value             = local_force->evaluate(point_entity);
        const auto diffusive_flux_gradient = local_diffusive_flux_reconstruction->jacobian(point_entity);
        // compute the residual estimator
        const auto diffusive_flux_divergence = diffusive_flux_gradient[0][0] + diffusive_flux_gradient[1][1];
        const auto residual_difference       = force_value - diffusive_flux_divergence;
        const auto residual_product          = residual_difference * residual_difference;
        estimators_residual[entity_index] += integration_factor * quadrature_weight * residual_product;
      } // do a volume quadrature
      // for the residual estimator
      estimators_residual[entity_index] *=
          (entity.geometry().volume() * poincare_constant_sqrt) / std::sqrt(diffusion_min);
    } // walk the grid for the third time

    // compute error components
    double residual_estimator = 0;
    for (size_t ii = 0; ii < estimators_residual.size(); ++ii) {
      residual_estimator += std::pow(estimators_residual[ii], 2);
    }
    return std::sqrt(residual_estimator);
  } // ... compute_residual_estimator_ESV10(...)

  double compute_diffusive_flux_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const size_t integration_order = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(
        grid_part, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(
        discretization.space(), *current_solution_vector_on_level_, "discrete solution");

    typedef RaviartThomasSpace::FemLocalfunctionsWrapper<GridPartType, 0, RangeFieldType, dimDomain> RTN_SpaceType;
    const RTN_SpaceType rtn_space(grid_part);
    typedef DiscreteFunction<RTN_SpaceType, VectorType> RTN_DiscreteFunctionType;
    VectorType rtn_vector(rtn_space.mapper().size());
    RTN_DiscreteFunctionType diffusive_flux_reconstruction(rtn_space, rtn_vector, "diffusive flux reconstruction");
    // reconstruct
    typedef ReconstructionOperator::DiffusiveFlux<GridPartType, typename TestCase::DiffusionType>
        ReconstructionOperatorType;
    const ReconstructionOperatorType reconstruction_operator(*grid_part, test_.diffusion());
    reconstruction_operator.apply(discrete_solution, diffusive_flux_reconstruction);

    // walk the grid for the third time
    std::vector<RangeFieldType> estimators_diffusive_flux(grid_part->indexSet().size(0), RangeFieldType(0));
    const auto entity_it_end = grid_part->template end<0>();
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion                     = test_.diffusion().local_function(entity);
      const auto local_solution                      = discrete_solution.local_function(entity);
      const auto local_diffusive_flux_reconstruction = diffusive_flux_reconstruction.local_function(entity);
      // do a volume quadrature
      const auto& volume_quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector<DomainFieldType, dimDomain> point_entity = quadrature_point.position();
        const double quadrature_weight                             = quadrature_point.weight();
        // evaluate
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value      = local_diffusion->evaluate(point_entity);
        auto solution_gradient          = local_solution->jacobian(point_entity);
        auto diffusive_flux_value       = local_diffusive_flux_reconstruction->evaluate(point_entity);
        // compute diffusive flux estimator
        const auto diffusion_value_sqrt = std::sqrt(diffusion_value);
        solution_gradient[0] *= diffusion_value_sqrt;
        diffusive_flux_value /= diffusion_value_sqrt;
        const auto diffusive_flux_sum     = solution_gradient[0] + diffusive_flux_value;
        const auto diffusive_flux_product = diffusive_flux_sum * diffusive_flux_sum;
        estimators_diffusive_flux[entity_index] += integration_factor * quadrature_weight * diffusive_flux_product;
      } // do a volume quadrature
    } // walk the grid for the third time

    // compute error components
    double diffusive_flux_estimator = 0;
    for (size_t ii = 0; ii < estimators_diffusive_flux.size(); ++ii) {
      diffusive_flux_estimator += std::pow(estimators_diffusive_flux[ii], 2);
    }
    return std::sqrt(diffusive_flux_estimator);
  } // ... compute_diffusive_flux_estimator_ESV10(...)

  double compute_nonconformity_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const size_t integration_order = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(
        grid_part, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(
        discretization.space(), *current_solution_vector_on_level_, "discrete solution");
    VectorType oswald_projection_vector(discretization.space().mapper().size());
    DiscreteFunctionType oswald_projection(discretization.space(), oswald_projection_vector, "oswald projection");

    typedef typename DiscretizationType::SpaceType TestSpaceType;
    typedef FieldVector<DomainFieldType, dimDomain> DomainType;

    // data structures we need
    // * a map from a global vertex id to global DoF ids
    //   given a vertex, one obtains a set of all global DoF ids, which are associated with this vertex
    typedef std::vector<std::set<size_t>> VertexToEntitiesMapType;
    VertexToEntitiesMapType vertex_to_dof_id_map(grid_part->indexSet().size(dimDomain));
    // * a set to hold the global id off all boundary vertices
    std::set<size_t> boundary_vertices;
    // * vectors to hold the local estimators
    std::vector<RangeFieldType> estimators_nonconformity(grid_part->indexSet().size(0), RangeFieldType(0));

    // walk the grid for the first time
    const auto entity_it_end = grid_part->template end<0>();
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local finite elements
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // loop over all vertices
      std::vector<DomainType> global_vertices(entity.template count<dimDomain>(), DomainType(0));
      std::vector<size_t> global_vertex_ids(global_vertices.size(), 0);
      assert(global_vertices.size() < std::numeric_limits<int>::max());
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get global vertex id
        const auto vertexPtr               = entity.template subEntity<dimDomain>(int(local_vertex_id));
        const auto& vertex                 = *vertexPtr;
        global_vertex_ids[local_vertex_id] = grid_part->indexSet().index(vertex);
        global_vertices[local_vertex_id]   = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            const size_t local_DOF_id  = entity_dg_local_key.index();
            const size_t global_DOF_id = discretization.space().mapper().mapToGlobal(entity, local_DOF_id);
            // add this global DoF to this vertex
            vertex_to_dof_id_map[global_vertex_ids[local_vertex_id]].insert(global_DOF_id);
            // there must be one and only one for a polorder 1 lagrange basis
            break;
          }
        } // loop over all local DoFs
      } // loop over all vertices
      // in order to determine the boundary vertices, we need to
      // loop over all intersections
      const auto intersectionEndIt = grid_part->iend(entity);
      for (auto intersectionIt = grid_part->ibegin(entity); intersectionIt != intersectionEndIt; ++intersectionIt) {
        const auto& intersection = *intersectionIt;
        if (intersection.boundary() && !intersection.neighbor()) {
          const auto& intersection_geometry = intersection.geometry();
          for (size_t local_intersection_corner_id = 0;
               int(local_intersection_corner_id) < intersection_geometry.corners();
               ++local_intersection_corner_id) {
            const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
            // now, we need to find the entity's vertex this intersection's corner point equals to, so we
            // loop over all vertices of the entity
            for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
              if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id]))
                boundary_vertices.insert(global_vertex_ids[local_vertex_id]);
            } // loop over all vertices of the entity
          } // if (intersection.boundary() && !intersection.neighbor())
        } // loop over all intersections
      } // loop over all intersections
    } // walk the grid for the first time

    // walk the grid for the second time
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_solution_entity             = discrete_solution.local_discrete_function(entity);
      const auto& local_solution_entity_DoF_vector = local_solution_entity.vector();
      // get the local finite elements
      // * for the oswald projection
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // to compute the oswald projection
      // * loop over all local DoFs
      for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
        const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
        const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
        assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
        const size_t local_vertex_id  = entity_cg_local_key.subEntity();
        const size_t local_DoF_id     = entity_dg_local_key.index();
        const auto vertexPtr          = entity.template subEntity<dimDomain>(local_vertex_id);
        const auto& vertex            = *vertexPtr;
        const size_t global_vertex_id = grid_part->indexSet().index(vertex);
        // if we are on the domain boundary
        if (boundary_vertices.count(global_vertex_id)) {
          // get global DoF id
          const size_t global_DoF_id = discretization.space().mapper().mapToGlobal(entity, local_DoF_id);
          // set the dof to zero (we have dirichlet zero)
          oswald_projection.vector().set(global_DoF_id, RangeFieldType(0));
        } else {
          // do the oswald projection
          const size_t num_DoFS_per_vertex = vertex_to_dof_id_map[global_vertex_id].size();
          // * get the source Dof
          const RangeFieldType source_Dof_value = local_solution_entity_DoF_vector.get(local_DoF_id);
          // * and add it to all target DoFs
          for (size_t target_global_DoF_id : vertex_to_dof_id_map[global_vertex_id])
            oswald_projection.vector().add(target_global_DoF_id, source_Dof_value / num_DoFS_per_vertex);
        } // if (boundary_vertices.find(global_vertex_id))
      } // loop over all local DoFs
    } // walk the grid for the second time

    // walk the grid for the third time
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion         = test_.diffusion().local_function(entity);
      const auto local_solution          = discrete_solution.local_function(entity);
      const auto local_oswald_projection = oswald_projection.local_function(entity);
      // do a volume quadrature
      const auto& volume_quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector<DomainFieldType, dimDomain> point_entity = quadrature_point.position();
        const double quadrature_weight                             = quadrature_point.weight();
        // evaluate
        const double integration_factor       = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value            = local_diffusion->evaluate(point_entity);
        auto solution_gradient                = local_solution->jacobian(point_entity);
        const auto oswald_projection_gradient = local_oswald_projection->jacobian(point_entity);
        // compute local nonconformity estimator
        const auto nonconformity_difference = solution_gradient[0] - oswald_projection_gradient[0];
        const auto nonconformity_product    = nonconformity_difference * nonconformity_difference;
        estimators_nonconformity[entity_index] +=
            integration_factor * quadrature_weight * diffusion_value * nonconformity_product;
      } // do a volume quadrature
    } // walk the grid for the third time

    // compute error components
    double nonconformity_estimator = 0;
    for (size_t ii = 0; ii < estimators_nonconformity.size(); ++ii) {
      nonconformity_estimator += std::pow(estimators_nonconformity[ii], 2);
    }
    return std::sqrt(nonconformity_estimator);
  }

  double compute_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const auto poincare_constant_sqrt = RangeFieldType(1) / M_PIl;
    const size_t integration_order    = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(
        grid_part, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(
        discretization.space(), *current_solution_vector_on_level_, "discrete solution");
    VectorType oswald_projection_vector(discretization.space().mapper().size());
    DiscreteFunctionType oswald_projection(discretization.space(), oswald_projection_vector, "oswald projection");

    typedef typename DiscretizationType::SpaceType TestSpaceType;
    typedef FieldVector<DomainFieldType, dimDomain> DomainType;

    typedef RaviartThomasSpace::FemLocalfunctionsWrapper<GridPartType, 0, RangeFieldType, dimDomain> RTN_SpaceType;
    const RTN_SpaceType rtn_space(grid_part);
    typedef DiscreteFunction<RTN_SpaceType, VectorType> RTN_DiscreteFunctionType;
    VectorType rtn_vector(rtn_space.mapper().size());
    RTN_DiscreteFunctionType diffusive_flux_reconstruction(rtn_space, rtn_vector, "diffusive flux reconstruction");
    // reconstruct
    typedef ReconstructionOperator::DiffusiveFlux<GridPartType, typename TestCase::DiffusionType>
        ReconstructionOperatorType;
    const ReconstructionOperatorType reconstruction_operator(*grid_part, test_.diffusion());
    reconstruction_operator.apply(discrete_solution, diffusive_flux_reconstruction);

    // data structures we need
    // * a map from a global vertex id to global DoF ids
    //   given a vertex, one obtains a set of all global DoF ids, which are associated with this vertex
    typedef std::vector<std::set<size_t>> VertexToEntitiesMapType;
    VertexToEntitiesMapType vertex_to_dof_id_map(grid_part->indexSet().size(dimDomain));
    // * a set to hold the global id off all boundary vertices
    std::set<size_t> boundary_vertices;
    // * vectors to hold the local estimators
    std::vector<RangeFieldType> estimators_nonconformity(grid_part->indexSet().size(0), RangeFieldType(0));
    std::vector<RangeFieldType> estimators_residual(grid_part->indexSet().size(0), RangeFieldType(0));
    std::vector<RangeFieldType> estimators_diffusive_flux(grid_part->indexSet().size(0), RangeFieldType(0));

    // walk the grid for the first time
    const auto entity_it_end = grid_part->template end<0>();
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local finite elements
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // loop over all vertices
      std::vector<DomainType> global_vertices(entity.template count<dimDomain>(), DomainType(0));
      std::vector<size_t> global_vertex_ids(global_vertices.size(), 0);
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get global vertex id
        const auto vertexPtr               = entity.template subEntity<dimDomain>(local_vertex_id);
        const auto& vertex                 = *vertexPtr;
        global_vertex_ids[local_vertex_id] = grid_part->indexSet().index(vertex);
        global_vertices[local_vertex_id]   = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            const size_t local_DOF_id  = entity_dg_local_key.index();
            const size_t global_DOF_id = discretization.space().mapper().mapToGlobal(entity, local_DOF_id);
            // add this global DoF to this vertex
            vertex_to_dof_id_map[global_vertex_ids[local_vertex_id]].insert(global_DOF_id);
            // there must be one and only one for a polorder 1 lagrange basis
            break;
          }
        } // loop over all local DoFs
      } // loop over all vertices
      // in order to determine the boundary vertices, we need to
      // loop over all intersections
      const auto intersectionEndIt = grid_part->iend(entity);
      for (auto intersectionIt = grid_part->ibegin(entity); intersectionIt != intersectionEndIt; ++intersectionIt) {
        const auto& intersection = *intersectionIt;
        if (intersection.boundary() && !intersection.neighbor()) {
          const auto& intersection_geometry = intersection.geometry();
          for (size_t local_intersection_corner_id = 0;
               int(local_intersection_corner_id) < intersection_geometry.corners();
               ++local_intersection_corner_id) {
            const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
            // now, we need to find the entity's vertex this intersection's corner point equals to, so we
            // loop over all vertices of the entity
            for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
              if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id]))
                boundary_vertices.insert(global_vertex_ids[local_vertex_id]);
            } // loop over all vertices of the entity
          } // if (intersection.boundary() && !intersection.neighbor())
        } // loop over all intersections
      } // loop over all intersections
    } // walk the grid for the first time

    // walk the grid for the second time
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_solution_entity             = discrete_solution.local_discrete_function(entity);
      const auto& local_solution_entity_DoF_vector = local_solution_entity.vector();
      // get the local finite elements
      // * for the oswald projection
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // * for the diffusive flux reconstruction
      const auto rtn_finite_element      = rtn_space.backend().finiteElement(entity);
      const auto& rtn_local_coefficients = rtn_finite_element.localCoefficients();
      std::vector<size_t> intersection_id_to_local_DoF_id_map(rtn_local_coefficients.size(), 0);
      for (size_t ii = 0; ii < rtn_local_coefficients.size(); ++ii)
        intersection_id_to_local_DoF_id_map[rtn_local_coefficients.localKey(ii).subEntity()] = ii;
      // to compute the oswald projection
      // * loop over all local DoFs
      for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
        const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
        const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
        assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
        const size_t local_vertex_id  = entity_cg_local_key.subEntity();
        const size_t local_DoF_id     = entity_dg_local_key.index();
        const auto vertexPtr          = entity.template subEntity<dimDomain>(local_vertex_id);
        const auto& vertex            = *vertexPtr;
        const size_t global_vertex_id = grid_part->indexSet().index(vertex);
        // if we are on the domain boundary
        if (boundary_vertices.count(global_vertex_id)) {
          // get global DoF id
          const size_t global_DoF_id = discretization.space().mapper().mapToGlobal(entity, local_DoF_id);
          // set the dof to zero (we have dirichlet zero)
          oswald_projection.vector().set(global_DoF_id, RangeFieldType(0));
        } else {
          // do the oswald projection
          const size_t num_DoFS_per_vertex = vertex_to_dof_id_map[global_vertex_id].size();
          // * get the source Dof
          const RangeFieldType source_Dof_value = local_solution_entity_DoF_vector.get(local_DoF_id);
          // * and add it to all target DoFs
          for (size_t target_global_DoF_id : vertex_to_dof_id_map[global_vertex_id])
            oswald_projection.vector().add(target_global_DoF_id, source_Dof_value / num_DoFS_per_vertex);
        } // if (boundary_vertices.find(global_vertex_id))
      } // loop over all local DoFs
    } // walk the grid for the second time

    // walk the grid for the third time
    for (auto entity_it = grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion                     = test_.diffusion().local_function(entity);
      const auto local_force                         = test_.force().local_function(entity);
      const auto local_solution                      = discrete_solution.local_function(entity);
      const auto local_oswald_projection             = oswald_projection.local_function(entity);
      const auto local_diffusive_flux_reconstruction = diffusive_flux_reconstruction.local_function(entity);
      RangeFieldType diffusion_min                   = std::numeric_limits<RangeFieldType>::max();
      // do a volume quadrature
      const auto& volume_quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), 2 * integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector<DomainFieldType, dimDomain> point_entity = quadrature_point.position();
        const double quadrature_weight                             = quadrature_point.weight();
        // evaluate
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value = local_diffusion->evaluate(point_entity);
        if (diffusion_value < diffusion_min)
          diffusion_min                       = diffusion_value;
        const auto force_value                = local_force->evaluate(point_entity);
        auto solution_gradient                = local_solution->jacobian(point_entity);
        const auto oswald_projection_gradient = local_oswald_projection->jacobian(point_entity);
        auto diffusive_flux_value             = local_diffusive_flux_reconstruction->evaluate(point_entity);
        const auto diffusive_flux_gradient    = local_diffusive_flux_reconstruction->jacobian(point_entity);
        // compute local nonconformity estimator
        const auto nonconformity_difference = solution_gradient[0] - oswald_projection_gradient[0];
        const auto nonconformity_product    = nonconformity_difference * nonconformity_difference;
        estimators_nonconformity[entity_index] +=
            integration_factor * quadrature_weight * diffusion_value * nonconformity_product;
        // compute the residual estimator
        const auto diffusive_flux_divergence = diffusive_flux_gradient[0][0] + diffusive_flux_gradient[1][1];
        const auto residual_difference       = force_value - diffusive_flux_divergence;
        const auto residual_product          = residual_difference * residual_difference;
        estimators_residual[entity_index] += integration_factor * quadrature_weight * residual_product;
        // compute diffusive flux estimator
        const auto diffusion_value_sqrt = std::sqrt(diffusion_value);
        solution_gradient[0] *= diffusion_value_sqrt;
        diffusive_flux_value /= diffusion_value_sqrt;
        const auto diffusive_flux_sum     = solution_gradient[0] + diffusive_flux_value;
        const auto diffusive_flux_product = diffusive_flux_sum * diffusive_flux_sum;
        estimators_diffusive_flux[entity_index] += integration_factor * quadrature_weight * diffusive_flux_product;
      } // do a volume quadrature
      // for the residual estimator
      estimators_residual[entity_index] *=
          (entity.geometry().volume() * poincare_constant_sqrt) / std::sqrt(diffusion_min);
    } // walk the grid for the third time

    // compute error components
    double estimator = 0;
    assert(estimators_nonconformity.size() == estimators_residual.size());
    assert(estimators_nonconformity.size() == estimators_diffusive_flux.size());
    for (size_t ii = 0; ii < estimators_nonconformity.size(); ++ii) {
      estimator += std::pow(estimators_nonconformity[ii], 2)
                   + std::pow(estimators_residual[ii] + estimators_diffusive_flux[ii], 2);
    }
    return std::sqrt(estimator);
  } // ... compute_estimator_ESV10(...)

private:
  using BaseType::test_;
  using BaseType::current_level_;
  using BaseType::last_computed_level_;
  using BaseType::reference_solution_computed_;
  using BaseType::reference_discretization_;
  using BaseType::reference_solution_vector_;
  using BaseType::current_solution_vector_on_level_;
  using BaseType::current_solution_vector_;
}; // class EstimatorStudy


} // namespace EllipticSWIPDG

#endif // DUNE_GDT_TEST_ELLIPTIC_SWIPDG_DISCRETIZATION_HH
