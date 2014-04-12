// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_SIPDG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_SIPDG_DISCRETIZATION_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/common/timer.hh>

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
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localevaluation/sipdg.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/space/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/product/l2.hh>
#include <dune/gdt/product/h1.hh>
#include <dune/gdt/operator/prolongations.hh>

#include "elliptic-testcases.hh"


namespace EllipticSIPDG {


template <class GridPartType, int polynomialOrder, class MatrixImp = Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,
          class VectorImp                                          = Dune::Stuff::LA::EigenDenseVector<double>>
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

  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;

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
    , is_assembled_(false)
    , system_matrix_(0, 0)
    , rhs_vector_(0)
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

  void assemble() const
  {
    if (!is_assembled_) {
      using namespace Dune;
      using namespace Dune::GDT;

      Stuff::LA::SparsityPatternDefault sparsity_pattern = space_.compute_face_and_volume_pattern();
      system_matrix_                                     = MatrixType(space_.mapper().size(), space_.mapper().size(), sparsity_pattern);
      rhs_vector_                                        = VectorType(space_.mapper().size());
      typedef SystemAssembler<SpaceType> SystemAssemblerType;
      SystemAssemblerType systemAssembler(space_);

      // volume terms
      // * lhs
      typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<FunctionType>> EllipticOperatorType;
      const EllipticOperatorType ellipticOperator(diffusion_);
      const LocalAssembler::Codim0Matrix<EllipticOperatorType> diffusionMatrixAssembler(ellipticOperator);
      systemAssembler.add(diffusionMatrixAssembler, system_matrix_);
      // * rhs
      typedef LocalFunctional::Codim0Integral<LocalEvaluation::Product<FunctionType>> ForceFunctionalType;
      const ForceFunctionalType forceFunctional(force_);
      const LocalAssembler::Codim0Vector<ForceFunctionalType> forceVectorAssembler(forceFunctional);
      systemAssembler.add(forceVectorAssembler, rhs_vector_);
      // inner face terms
      typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SIPDG::Inner<FunctionType>> CouplingOperatorType;
      const CouplingOperatorType couplingOperator(diffusion_, beta_);
      const LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> couplingMatrixAssembler(couplingOperator);
      systemAssembler.add(couplingMatrixAssembler,
                          system_matrix_,
                          new ApplyOn::InnerIntersectionsPrimally<typename SpaceType::GridViewType>());
      // dirichlet boundary face terms
      // * lhs
      typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SIPDG::BoundaryLHS<FunctionType>>
          DirichletOperatorType;
      const DirichletOperatorType dirichletOperator(diffusion_, beta_);
      const LocalAssembler::Codim1BoundaryMatrix<DirichletOperatorType> dirichletMatrixAssembler(dirichletOperator);
      systemAssembler.add(dirichletMatrixAssembler,
                          system_matrix_,
                          new ApplyOn::DirichletIntersections<typename SpaceType::GridViewType>(boundary_info_));
      // * rhs
      typedef LocalFunctional::Codim1Integral<LocalEvaluation::SIPDG::BoundaryRHS<FunctionType, FunctionType>>
          DirichletFunctionalType;
      const DirichletFunctionalType dirichletFunctional(diffusion_, dirichlet_, beta_);
      const LocalAssembler::Codim1Vector<DirichletFunctionalType> dirichletVectorAssembler(dirichletFunctional);
      systemAssembler.add(dirichletVectorAssembler,
                          rhs_vector_,
                          new ApplyOn::DirichletIntersections<typename SpaceType::GridViewType>(boundary_info_));
      // neumann boundary face terms
      // * rhs
      typedef LocalFunctional::Codim1Integral<LocalEvaluation::Product<FunctionType>> NeumannFunctionalType;
      const NeumannFunctionalType neumannFunctional(neumann_);
      const LocalAssembler::Codim1Vector<NeumannFunctionalType> neumannVectorAssembler(neumannFunctional);
      systemAssembler.add(neumannVectorAssembler,
                          rhs_vector_,
                          new ApplyOn::NeumannIntersections<typename SpaceType::GridViewType>(boundary_info_));
      // do all the work
      systemAssembler.assemble();
      is_assembled_ = true;
    }
  } // ... assemble()

  bool assembled() const
  {
    return is_assembled_;
  }

  const MatrixType& system_matrix() const
  {
    return system_matrix_;
  }

  const VectorType& rhs_vector() const
  {
    return rhs_vector_;
  }

  void solve(VectorType& solution) const
  {
    if (!is_assembled_)
      assemble();
    Dune::Stuff::LA::Solver<MatrixType>(system_matrix_).apply(rhs_vector_, solution);
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
  const RangeFieldType beta_;
  mutable bool is_assembled_;
  mutable MatrixType system_matrix_;
  mutable VectorType rhs_vector_;
}; // class Discretization


template <class TestCase, int polOrder>
class EocStudy : public Dune::Stuff::Common::ConvergenceStudy
{
  typedef Dune::Stuff::Common::ConvergenceStudy BaseType;
  typedef typename TestCase::GridPartType GridPartType;
  typedef typename TestCase::GridViewType GridViewType;
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

  virtual std::string identifier() const DS_OVERRIDE DS_FINAL
  {
    return "SIP discontinuous galerkin discretization, polOrder " + Dune::Stuff::Common::toString(polOrder);
  }

  virtual size_t num_refinements() const DS_OVERRIDE DS_FINAL
  {
    if (test_.num_levels() == 0)
      return test_.num_levels();
    else
      return test_.num_levels() - 1;
  }

  virtual std::vector<std::string> provided_norms() const DS_OVERRIDE DS_FINAL
  {
    return {"L2", "H1_semi"};
  }

  virtual size_t expected_rate(const std::string type) const DS_OVERRIDE DS_FINAL
  {
    if (type.compare("L2") == 0)
      return polOrder + 1;
    else if (type.compare("H1_semi") == 0)
      return polOrder;
    else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  virtual double norm_reference_solution(const std::string type) DS_OVERRIDE DS_FINAL
  {
    if (test_.provides_exact_solution()) {
      return compute_norm(*(test_.reference_grid_view()), test_.exact_solution(), type);
    } else {
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "reference solution");
      // compute norm
      return compute_norm(*(test_.reference_grid_view()), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const DS_OVERRIDE DS_FINAL
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return grid_part->grid().size(grid_part->level(), 0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const DS_OVERRIDE DS_FINAL
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_part = test_.level_grid_part(current_level_);
    return Dune::Fem::GridWidth::calcGridWidth(*grid_part);
  } // ... current_grid_width(...)

  virtual double compute_on_current_refinement() DS_OVERRIDE DS_FINAL
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
      auto current_level_solution_vector = discretization.create_vector();
      discretization.solve(current_level_solution_vector);
      const ConstDiscreteFunctionType current_level_solution(
          discretization.space(), current_level_solution_vector, "solution on current level");
      // prolong to reference grid part
      elapsed += timer.elapsed();
      if (!reference_solution_computed_)
        compute_reference_solution();
      timer.reset();
      const auto reference_grid_view = test_.reference_grid_view();
      assert(reference_discretization_);
      if (!current_solution_vector_)
        current_solution_vector_ =
            std::unique_ptr<VectorType>(new VectorType(reference_discretization_->create_vector()));
      DiscreteFunctionType reference_level_solution(
          reference_discretization_->space(), *current_solution_vector_, "solution on reference grid view");
      const ProlongationOperator::L2<GridViewType> prolongation_operator(*reference_grid_view);
      prolongation_operator.apply(current_level_solution, reference_level_solution);
      last_computed_level_ = current_level_;
    }
    return timer.elapsed() + elapsed;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) DS_OVERRIDE DS_FINAL
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
      return compute_norm(*(test_.reference_grid_view()), difference, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "reference solution");
      typedef Dune::Stuff::Function::Difference<ConstDiscreteFunctionType, ConstDiscreteFunctionType> DifferenceType;
      const DifferenceType difference(reference_solution, current_solution);
      return compute_norm(*(test_.reference_grid_view()), difference, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() DS_OVERRIDE DS_FINAL
  {
    if (current_level_ < test_.num_levels())
      ++current_level_;
  }

  std::vector<double> expected_results(const std::string type) const
  {
    if (std::is_same<TestCase, EllipticTestCase::ESV07<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.24e-01, 3.43e-02, 8.83e-03, 2.24e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.73e-01, 1.91e-01, 9.55e-02, 4.79e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {1.38e-02, 1.75e-03, 2.19e-04, 2.74e-05};
        else if (type.compare("H1_semi") == 0)
          return {7.43e-02, 1.92e-02, 4.84e-03, 1.22e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::LocalThermalBlock<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {6.80e-02, 3.87e-02, 1.11e-02, 2.47e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.48e-01, 3.00e-01, 1.63e-01, 7.47e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {2.48e-02, 3.72e-03, 6.29e-04, 1.22e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.09e-01, 6.94e-02, 2.21e-02, 6.91e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::ER07<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {8.43e-02, 2.26e-02, 5.75e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.07e-01, 1.56e-01, 7.75e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {7.28e-03, 9.41e-04, 1.20e-04};
        else if (type.compare("H1_semi") == 0)
          return {4.97e-02, 1.28e-02, 3.21e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::MixedBoundaryTypes<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {4.76e-02, 1.36e-02, 3.50e-03, 7.77e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.74e-01, 1.45e-01, 7.16e-02, 3.23e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {3.96e-03, 7.11e-04, 1.44e-04, 3.24e-05};
        else if (type.compare("H1_semi") == 0)
          return {4.90e-02, 1.83e-02, 7.39e-03, 2.97e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  }

  std::map<std::string, std::vector<double>> run(std::ostream& out = std::cout)
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
      if (!reference_solution_vector_)
        reference_solution_vector_ =
            std::unique_ptr<VectorType>(new VectorType(reference_discretization_->create_vector()));
      reference_discretization_->solve(*reference_solution_vector_);
      reference_solution_computed_ = true;
    }
  } // ... compute_reference_solution()

  template <class GridViewType, class FunctionType>
  double compute_norm(const GridViewType& grid_view, const FunctionType& function, const std::string type)
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (type.compare("L2") == 0) {
      Product::L2<GridViewType> l2_product_operator(grid_view);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      Product::H1SemiGeneric<GridViewType> h1_product_operator(grid_view);
      return std::sqrt(h1_product_operator.apply2(function, function));
    } else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... compute_norm(...)

  const TestCase& test_;
  size_t current_level_;
  size_t last_computed_level_;
  bool reference_solution_computed_;
  std::unique_ptr<DiscretizationType> reference_discretization_;
  std::unique_ptr<VectorType> reference_solution_vector_;
  std::unique_ptr<VectorType> current_solution_vector_;
}; // class EOCStudy


} // namespace EllipticSIPDG

#endif // DUNE_GDT_TEST_ELLIPTIC_SIPDG_DISCRETIZATION_HH
