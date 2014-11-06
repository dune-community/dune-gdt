// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_SIPDG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_SIPDG_DISCRETIZATION_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/common/timer.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/playground/spaces/discontinuouslagrange/fem.hh>
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
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/operators/prolongations.hh>

#include "elliptic-testcases.hh"

#if HAVE_EIGEN
static auto constexpr matrix_tag = Dune::Stuff::LA::ChooseBackend::eigen_sparse;
#else
static auto constexpr matrix_tag = Dune::Stuff::LA::default_sparse_backend;
#endif

namespace EllipticSIPDG {


template <class GridPartType, int polynomialOrder,
          class MatrixImp = typename Dune::Stuff::LA::Container<double, matrix_tag>::MatrixType,
          class VectorImp = typename Dune::Stuff::LA::Container<double, matrix_tag>::VectorType>
class Discretization
{
public:
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype DomainFieldType;

  static const unsigned int dimRange = 1;
  typedef double RangeFieldType;

  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::Grid::BoundaryInfoInterface<typename GridPartType::IntersectionType> BoundaryInfoType;
  typedef Dune::Stuff::LocalizableFunctionInterface<typename GridPartType::template Codim<0>::EntityType,
                                                    DomainFieldType, dimDomain, RangeFieldType, dimRange> FunctionType;

  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;

  typedef Dune::GDT::Spaces::DiscontinuousLagrange::FemBased<GridPartType, polOrder, RangeFieldType, dimRange>
      SpaceType;

  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;

  Discretization(const GridPartType& gp, const BoundaryInfoType& info, const FunctionType& diff,
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
                          new Stuff::Grid::ApplyOn::InnerIntersectionsPrimally<typename SpaceType::GridViewType>());
      // dirichlet boundary face terms
      // * lhs
      typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SIPDG::BoundaryLHS<FunctionType>>
          DirichletOperatorType;
      const DirichletOperatorType dirichletOperator(diffusion_, beta_);
      const LocalAssembler::Codim1BoundaryMatrix<DirichletOperatorType> dirichletMatrixAssembler(dirichletOperator);
      systemAssembler.add(
          dirichletMatrixAssembler,
          system_matrix_,
          new Stuff::Grid::ApplyOn::DirichletIntersections<typename SpaceType::GridViewType>(boundary_info_));
      // * rhs
      typedef LocalFunctional::Codim1Integral<LocalEvaluation::SIPDG::BoundaryRHS<FunctionType, FunctionType>>
          DirichletFunctionalType;
      const DirichletFunctionalType dirichletFunctional(diffusion_, dirichlet_, beta_);
      const LocalAssembler::Codim1Vector<DirichletFunctionalType> dirichletVectorAssembler(dirichletFunctional);
      systemAssembler.add(
          dirichletVectorAssembler,
          rhs_vector_,
          new Stuff::Grid::ApplyOn::DirichletIntersections<typename SpaceType::GridViewType>(boundary_info_));
      // neumann boundary face terms
      // * rhs
      typedef LocalFunctional::Codim1Integral<LocalEvaluation::Product<FunctionType>> NeumannFunctionalType;
      const NeumannFunctionalType neumannFunctional(neumann_);
      const LocalAssembler::Codim1Vector<NeumannFunctionalType> neumannVectorAssembler(neumannFunctional);
      systemAssembler.add(
          neumannVectorAssembler,
          rhs_vector_,
          new Stuff::Grid::ApplyOn::NeumannIntersections<typename SpaceType::GridViewType>(boundary_info_));
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

  virtual std::string identifier() const override final
  {
    return "SIP discontinuous galerkin discretization, polOrder " + Dune::Stuff::Common::toString(polOrder);
  }

  virtual size_t num_refinements() const override final
  {
    if (test_.num_levels() == 0)
      return test_.num_levels();
    else
      return test_.num_levels() - 1;
  }

  virtual std::vector<std::string> provided_norms() const override final
  {
    return {"L2", "H1_semi"};
  }

  virtual size_t expected_rate(const std::string type) const override final
  {
    if (type.compare("L2") == 0)
      return polOrder + 1;
    else if (type.compare("H1_semi") == 0)
      return polOrder;
    else
      DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
  } // ... expected_rate(...)

  virtual double norm_reference_solution(const std::string type) override final
  {
    if (test_.provides_exact_solution()) {
      return compute_norm(test_.reference_grid_view(), test_.exact_solution(), type);
    } else {
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "reference solution");
      // compute norm
      return compute_norm(test_.reference_grid_view(), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const override final
  {
    assert(current_level_ < test_.num_levels());
    return test_.level_grid_part(current_level_).indexSet().size(0);
  }

  virtual double compute_on_current_refinement() override final
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
      const Operators::L2Prolongation<GridViewType> prolongation_operator(reference_grid_view);
      prolongation_operator.apply(current_level_solution, reference_level_solution);
      last_computed_level_ = current_level_;
    }
    return timer.elapsed() + elapsed;
  } // ... compute_on_current_refinement(...)

  virtual double current_error_norm(const std::string type) override final
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
      typedef Dune::Stuff::Functions::Difference<ExactSolutionType, ConstDiscreteFunctionType> DifferenceType;
      const DifferenceType difference(test_.exact_solution(), current_solution);
      return compute_norm(test_.reference_grid_view(), difference, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "reference solution");
      typedef Dune::Stuff::Functions::Difference<ConstDiscreteFunctionType, ConstDiscreteFunctionType> DifferenceType;
      const DifferenceType difference(reference_solution, current_solution);
      return compute_norm(test_.reference_grid_view(), difference, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override final
  {
    if (current_level_ < test_.num_levels())
      ++current_level_;
  }

  std::vector<double> expected_results(const std::string type) const
  {
#if HAVE_ALUGRID
    if (std::is_same<TestCase, EllipticTestCase::ESV07<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>::value) {
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
    } else if (std::is_same<TestCase,
                            EllipticTestCase::LocalThermalBlock<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>::
                   value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {6.79e-02, 3.86e-02, 1.10e-02, 2.46e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.47e-01, 2.99e-01, 1.62e-01, 7.46e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {2.47e-02, 3.71e-03, 6.28e-04, 1.21e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.08e-01, 6.93e-02, 2.20e-02, 6.90e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase,
                            EllipticTestCase::ER07<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {8.42e-02, 2.25e-02, 5.74e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.06e-01, 1.55e-01, 7.74e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {7.28e-03, 9.40e-04, 1.19e-04};
        else if (type.compare("H1_semi") == 0)
          return {4.95e-02, 1.27e-02, 3.20e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase,
                            EllipticTestCase::
                                MixedBoundaryTypes<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {4.75e-02, 1.35e-02, 3.49e-03, 7.76e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.73e-01, 1.44e-01, 7.15e-02, 3.22e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {3.95e-03, 7.10e-04, 1.43e-04, 3.24e-05};
        else if (type.compare("H1_semi") == 0)
          return {4.89e-02, 1.82e-02, 7.38e-03, 2.96e-03};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
#endif
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
      Products::L2<GridViewType> l2_product_operator(grid_view);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      Products::H1Semi<GridViewType> h1_product_operator(grid_view);
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
