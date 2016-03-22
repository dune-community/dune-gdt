// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_SWIPDG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_SWIPDG_DISCRETIZATION_HH

#include <memory>
#include <vector>
#include <type_traits>
#include <cmath>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/timer.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/ESV2007.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localevaluation/swipdg.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/operators/fluxreconstruction.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/products/ESV2007.hh>
#include <dune/gdt/playground/spaces/dg/fem.hh>
#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/weightedl2.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/spaces/cg/fem.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/spaces/rt/pdelab.hh>

#include "elliptic-testcases.hh"


#if HAVE_EIGEN
static auto constexpr la_backend = Dune::Stuff::LA::ChooseBackend::eigen_sparse;
#else
static auto constexpr la_backend = Dune::Stuff::LA::default_sparse_backend;
#endif


namespace EllipticSWIPDG {


template< class GridPartType,
          int polynomialOrder,
          class MatrixImp = typename Dune::Stuff::LA::Container< double, la_backend >::MatrixType,
          class VectorImp = typename Dune::Stuff::LA::Container< double, la_backend >::VectorType >
class Discretization
{
public:
  static const size_t                   dimDomain = GridPartType::dimension;
  typedef typename GridPartType::ctype  DomainFieldType;

  static const size_t       dimRange = 1;
  typedef double            RangeFieldType;

  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::Grid::BoundaryInfoInterface< typename GridPartType::IntersectionType > BoundaryInfoType;
  typedef Dune::Stuff::LocalizableFunctionInterface
      < typename GridPartType::template Codim< 0 >::EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange >
    FunctionType;

  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;

  typedef Dune::GDT::Spaces::DG::FemBased
      < GridPartType, polOrder, RangeFieldType, dimRange > SpaceType;

  typedef typename SpaceType::GridViewType GridViewType;

  typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType >      DiscreteFunctionType;
  typedef Dune::GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;

  Discretization(const GridPartType& gp,
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
    , beta_(1.0)
    , is_assembled_(false)
    , system_matrix_(0, 0)
    , rhs_vector_(0)
  {}

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

      system_matrix_ = MatrixType(space_.mapper().size(),
                                  space_.mapper().size(),
                                  space_.compute_face_and_volume_pattern());
      rhs_vector_ = VectorType (space_.mapper().size());
      typedef SystemAssembler< SpaceType > SystemAssemblerType;
      SystemAssemblerType systemAssembler(space_);

      // volume terms
      // * lhs
      typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< FunctionType > > EllipticOperatorType;
      const EllipticOperatorType                                  ellipticOperator(diffusion_);
      const LocalAssembler::Codim0Matrix< EllipticOperatorType >  diffusionMatrixAssembler(ellipticOperator);
      systemAssembler.add(diffusionMatrixAssembler, system_matrix_);
      // * rhs
      typedef LocalFunctional::Codim0Integral< LocalEvaluation::Product< FunctionType > > ForceFunctionalType;
      const ForceFunctionalType                                 forceFunctional(force_);
      const LocalAssembler::Codim0Vector< ForceFunctionalType > forceVectorAssembler(forceFunctional);
      systemAssembler.add(forceVectorAssembler, rhs_vector_);
      // inner face terms
      typedef LocalOperator::Codim1CouplingIntegral< LocalEvaluation::SWIPDG::Inner< FunctionType > > CouplingOperatorType;
      const CouplingOperatorType                                          couplingOperator(diffusion_, beta_);
      const LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType >  couplingMatrixAssembler(couplingOperator);
      systemAssembler.add(couplingMatrixAssembler,
                          system_matrix_,
                          new Stuff::Grid::ApplyOn::InnerIntersectionsPrimally< GridViewType >());
      // dirichlet boundary face terms
      // * lhs
      typedef LocalOperator::Codim1BoundaryIntegral< LocalEvaluation::SWIPDG::BoundaryLHS< FunctionType > >
          DirichletOperatorType;
      const DirichletOperatorType                                         dirichletOperator(diffusion_, beta_);
      const LocalAssembler::Codim1BoundaryMatrix< DirichletOperatorType > dirichletMatrixAssembler(dirichletOperator);
      systemAssembler.add(dirichletMatrixAssembler,
                          system_matrix_,
                          new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(boundary_info_));
      // * rhs
      typedef LocalFunctional::Codim1Integral< LocalEvaluation::SWIPDG::BoundaryRHS< FunctionType, FunctionType > >
          DirichletFunctionalType;
      const DirichletFunctionalType                                 dirichletFunctional(diffusion_,
                                                                                        dirichlet_,
                                                                                        beta_);
      const LocalAssembler::Codim1Vector< DirichletFunctionalType > dirichletVectorAssembler(dirichletFunctional);
      systemAssembler.add(dirichletVectorAssembler,
                          rhs_vector_,
                          new Stuff::Grid::ApplyOn::DirichletIntersections< GridViewType >(boundary_info_));
      // neumann boundary face terms
      // * rhs
      typedef LocalFunctional::Codim1Integral< LocalEvaluation::Product< FunctionType > > NeumannFunctionalType;
      const NeumannFunctionalType                                 neumannFunctional(neumann_);
      const LocalAssembler::Codim1Vector< NeumannFunctionalType > neumannVectorAssembler(neumannFunctional);
      systemAssembler.add(neumannVectorAssembler,
                          rhs_vector_,
                          new Stuff::Grid::ApplyOn::NeumannIntersections< GridViewType >(boundary_info_));
      // do all the work
      systemAssembler.assemble();
      is_assembled_ = true;
    }
  } // ... solve()

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
    Dune::Stuff::LA::Solver< MatrixType >(system_matrix_).apply(rhs_vector_, solution);
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
  mutable bool is_assembled_;
  mutable MatrixType system_matrix_;
  mutable VectorType rhs_vector_;
}; // class Discretization


template< class TestCase, int polOrder >
class EocStudy
  : public Dune::Stuff::Common::ConvergenceStudy
{
protected:
  typedef Dune::Stuff::Common::ConvergenceStudy BaseType;

  typedef typename TestCase::GridPartType GridPartType;
  typedef typename TestCase::EntityType   EntityType;
  typedef typename TestCase::GridViewType GridViewType;

  typedef typename TestCase::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = TestCase::dimDomain;
  typedef typename TestCase::RangeFieldType RangeFieldType;
  static const size_t                       dimRange = TestCase::dimRange;

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

  virtual std::string identifier() const override
  {
    return "SWIP discontinuous galerkin discretization, polOrder " + Dune::Stuff::Common::to_string(polOrder);
  }

  virtual size_t num_refinements() const override
  {
    if (test_.num_levels() == 0)
      return test_.num_levels();
    else
      return test_.num_levels() - 1;
  }

  virtual std::vector< std::string > provided_norms() const override
  {
    return {"L2", "H1_semi"};
  }

  virtual size_t expected_rate(const std::string type) const override
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
      return compute_norm(test_.reference_grid_view(), test_.exact_solution(), type);
    } else {
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(reference_discretization_->space(),
                                                         *reference_solution_vector_,
                                                         "reference solution");
      // compute norm
      return compute_norm(test_.reference_grid_view(), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const override
  {
    assert(current_level_ < test_.num_levels());
    return test_.level_grid_view(current_level_).indexSet().size(0);
  }

  virtual double current_grid_width() const override
  {
    assert(current_level_ < test_.num_levels());
    const GridViewType grid_view = test_.level_grid_view(current_level_);
    Dune::Stuff::Grid::Dimensions< GridViewType > dimensions(grid_view);
    return dimensions.entity_width.max();
  } // ... current_grid_width()

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
      const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                              test_.dirichlet(), test_.neumann());
      current_solution_vector_on_level_
          = std::unique_ptr< VectorType >(new VectorType(discretization.create_vector()));
      discretization.solve(*current_solution_vector_on_level_);
      const ConstDiscreteFunctionType current_level_solution(discretization.space(),
                                                             *current_solution_vector_on_level_,
                                                             "solution on current level");
      // prolong to reference grid part
      elapsed += timer.elapsed();
      if (!reference_solution_computed_)
        compute_reference_solution();
      timer.reset();
      const auto reference_grid_view = test_.reference_grid_view();
      const Operators::Prolongation< GridViewType > prolongation_operator(reference_grid_view);
      assert(reference_discretization_);
      if (!current_solution_vector_)
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
    const ConstDiscreteFunctionType current_solution(reference_discretization_->space(),
                                                     *current_solution_vector_,
                                                     "current solution");
    // compute error
    if (test_.provides_exact_solution()) {
      typedef Dune::Stuff::Functions::Difference< ExactSolutionType, ConstDiscreteFunctionType > DifferenceType;
      const DifferenceType difference(test_.exact_solution(), current_solution);
      return compute_norm(test_.reference_grid_view(), difference, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(reference_discretization_->space(),
                                                         *reference_solution_vector_,
                                                         "reference solution");
      typedef Dune::Stuff::Functions::Difference< ConstDiscreteFunctionType, ConstDiscreteFunctionType > DifferenceType;
      const DifferenceType difference(reference_solution, current_solution);
      return compute_norm(test_.reference_grid_view(), difference, type);
    }
  } // ... current_error_norm(...)

  virtual void refine() override
  {
    if (current_level_ < test_.num_levels())
      ++current_level_;
  }

  std::vector< double > expected_results(const std::string type) const override
  {
    using namespace Dune;
    if (std::is_same< TestCase, EllipticTestCase::ESV07< SGrid< 2, 2 > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.15e-01, 3.04e-02, 7.51e-03, 1.86e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.79e-01, 1.90e-01, 9.38e-02, 4.67e-02};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
#if HAVE_ALUGRID
    } else if (std::is_same< TestCase, EllipticTestCase::ESV07< ALUGrid<2, 2, Dune::simplex, Dune::conforming > > >::value
               || std::is_same< TestCase, EllipticTestCase::ESV07< ALUGrid< 2, 2, simplex, conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.82e-02, 4.53e-03, 1.12e-03, 2.78e-04};
        else if (type.compare("H1_semi") == 0)
          return {1.48e-01, 7.28e-02, 3.62e-02, 1.80e-02};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {8.55e-04, 1.06e-04, 1.31e-05, 1.63e-06};
        else if (type.compare("H1_semi") == 0)
          return {1.41e-02, 3.56e-03, 8.91e-04, 2.20e-04};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase,
               EllipticTestCase::LocalThermalBlock< ALUGrid<2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
#if THIS_IS_A_BUILDBOT_BUILD
        if (type.compare("L2") == 0)
          return {5.33e-02, 1.69e-02};
        else if (type.compare("H1_semi") == 0)
          return {3.82e-01, 2.29e-01};
#else // THIS_IS_A_BUILDBOT_BUILD
        if (type.compare("L2") == 0)
          return {5.57e-02, 1.99e-02, 5.54e-03, 1.29e-03};
        else if (type.compare("H1_semi") == 0)
          return {4.32e-01, 2.93e-01, 1.50e-01, 6.54e-02};
#endif // THIS_IS_A_BUILDBOT_BUILD
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
#if THIS_IS_A_BUILDBOT_BUILD
        if (type.compare("L2") == 0)
          return {1.18e-02, 2.12e-03};
        else if (type.compare("H1_semi") == 0)
          return {1.67e-01, 5.58e-02};
#else // THIS_IS_A_BUILDBOT_BUILD
        if (type.compare("L2") == 0)
          return {1.18e-02, 2.11e-03, 3.89e-04, 7.76e-05};
        else if (type.compare("H1_semi") == 0)
          return {1.69e-01, 5.96e-02, 1.94e-02, 6.04e-03};
#endif // THIS_IS_A_BUILDBOT_BUILD
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase,
               EllipticTestCase::ER07< ALUGrid<2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {6.09e-02, 1.65e-02, 4.22e-03};
        else if (type.compare("H1_semi") == 0)
          return {2.98e-01, 1.46e-01, 7.25e-02};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {6.42e-03, 8.23e-04, 1.04e-04};
        else if (type.compare("H1_semi") == 0)
          return {5.40e-02, 1.41e-02, 3.55e-03};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase,
               EllipticTestCase::MixedBoundaryTypes< ALUGrid<2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {4.02e-02, 1.12e-02, 2.83e-03, 6.33e-04};
        else if (type.compare("H1_semi") == 0)
          return {2.69e-01, 1.39e-01, 6.87e-02, 3.08e-02};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {3.58e-03, 6.25e-04, 1.21e-04, 2.68e-05};
        else if (type.compare("H1_semi") == 0)
          return {4.81e-02, 1.79e-02, 7.19e-03, 2.85e-03};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase,
               EllipticTestCase::Spe10Model1< ALUGrid<2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {7.22e-02, 2.59e-02};
        else if (type.compare("H1_semi") == 0)
          return {5.28e-01, 3.48e-01};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else if (polOrder == 2) {
        if (type.compare("L2") == 0)
          return {2.08e-02, 3.74e-03};
        else if (type.compare("H1_semi") == 0)
          return {2.56e-01, 8.47e-02};
        else
          DUNE_THROW(RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(NotImplemented, "Please record the expected results for this polOrder!");
#endif // HAVE_ALUGRID
    } else
      DUNE_THROW(NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  }

  std::map< std::string, std::vector< double > > run(std::ostream& out = std::cout)
  {
    return BaseType::run(true, out);
  }

protected:
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

  template< class GridViewType, class FunctionType >
  double compute_norm(const GridViewType& grid_view, const FunctionType& function, const std::string type)
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (type.compare("L2") == 0) {
      Products::L2< GridViewType > l2_product_operator(grid_view);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      Products::H1Semi< GridViewType > h1_product_operator(grid_view);
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
  std::unique_ptr< VectorType > current_solution_vector_on_level_;
  std::unique_ptr< VectorType > current_solution_vector_;
}; // class EOCStudy


template< class TestCase >
class EstimatorStudy
  : public EocStudy< TestCase, 1 >
{
  typedef EocStudy< TestCase, 1 > BaseType;
  static const unsigned int polOrder = 1;

  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType   EntityType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const size_t                         dimDomain = BaseType::dimDomain;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const size_t                       dimRange = BaseType::dimRange;

  typedef typename BaseType::DiscretizationType DiscretizationType;
  typedef typename DiscretizationType::VectorType                 VectorType;
  typedef typename DiscretizationType::DiscreteFunctionType       DiscreteFunctionType;
  typedef typename DiscretizationType::ConstDiscreteFunctionType  ConstDiscreteFunctionType;

  typedef typename TestCase::ExactSolutionType ExactSolutionType;

  static std::string nonconformity_estimator_id() {  return "eta_NC"; }
  static std::string residual_estimator_ESV07_id() { return "eta_R (ESV07)"; }
  static std::string diffusive_flux_estimator_id() { return "eta_DF"; }
  static std::string estimator_ESV07_id() {          return "eta (ESV07)"; }
  static std::string efficiency_ESV07_id() {         return "efficiency (ESV07)"; }

  const size_t over_integrate = 1;

public:
  EstimatorStudy(const TestCase& test)
    : BaseType(test)
  {}

  virtual ~EstimatorStudy() {}

  virtual std::vector< std::string > provided_norms() const override final
  {
    return { "energy"
            , nonconformity_estimator_id()
            , residual_estimator_ESV07_id()
            , diffusive_flux_estimator_id()
//            , estimator_ESV07_id()
            , efficiency_ESV07_id()
           };
  } // ... provided_norms(...)

  virtual size_t expected_rate(const std::string type) const override final
  {
    if (type == "energy")
      return polOrder;
    else if (type == nonconformity_estimator_id())
      return polOrder;
    else if (type == residual_estimator_ESV07_id())
      return polOrder + 1;
    else if (type == diffusive_flux_estimator_id())
      return polOrder;
    else if (type == estimator_ESV07_id())
      return polOrder;
    else if (type == efficiency_ESV07_id())
      return 0;
    else
      return BaseType::expected_rate(type);
  } // ... expected_rate(...)

  virtual double current_error_norm(const std::string type) override final
  {
    if (type == "energy")
      return compute_energy_norm();
    else if (type == nonconformity_estimator_id())
      return compute_nonconformity_estimator();
    else if (type == residual_estimator_ESV07_id())
      return compute_residual_estimator_ESV07();
    else if (type == diffusive_flux_estimator_id())
      return compute_diffusive_flux_estimator();
    else if (type == estimator_ESV07_id())
      return compute_estimator_ESV07();
    else if (type == efficiency_ESV07_id())
      return compute_efficiency_ESV07();
    else
      return BaseType::current_error_norm(type);
  } // ... current_error_norm(...)

  std::vector< double > expected_results(const std::string type) const override
  {
#if HAVE_ALUGRID
    if (std::is_same< TestCase, EllipticTestCase::ESV07< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("energy") == 0)
          return {3.28e-01, 1.62e-01, 8.04e-02, 4.01e-02};
        else if (type == nonconformity_estimator_id())
          return {1.66e-01, 7.89e-02, 3.91e-02, 1.95e-02};
        else if (type == residual_estimator_ESV07_id())
          return {7.23e-02, 1.82e-02, 4.54e-03, 1.14e-03};
        else if (type == diffusive_flux_estimator_id())
          return {3.55e-01, 1.76e-01, 8.73e-02, 4.35e-02};
//          return {3.39e-1, 1.70e-1, 8.40e-2, 4.19e-2};
        else if (type == efficiency_ESV07_id())
          return {1.37, 1.28, 1.23, 1.21};
//          return {1.21, 1.21, 1.21, 1.21};
        else
          return BaseType::expected_results(type);
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::LocalThermalBlock< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("energy") == 0)
#if THIS_IS_A_BUILDBOT_BUILD
          return {8.27e-02, 4.09e-02};
#else
          return {9.10e-02, 5.23e-02, 2.68e-02, 1.20e-02};
#endif
        else if (type == nonconformity_estimator_id())
          return {9.57e-02, 1.10e-01, 5.12e-02, 2.17e-02};
        else if (type == residual_estimator_ESV07_id())
          return {0.0, 0.0, 0.0, 0.0};
        else if (type == diffusive_flux_estimator_id())
          return {1.12e-01, 6.54e-02, 3.54e-02, 1.90e-02};
        else if (type == efficiency_ESV07_id())
#if THIS_IS_A_BUILDBOT_BUILD
          return {1.78e+00, 3.13e+00};
#else
          return {1.62e+00, 2.45e+00, 2.32e+00, 2.40e+00};
#endif
        else
          return BaseType::expected_results(type);
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same< TestCase, EllipticTestCase::Spe10Model1< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > > >::value) {
      if (polOrder == 1) {
        if (type.compare("energy") == 0)
          return {3.30e-02, 1.65e-02};
        else if (type == nonconformity_estimator_id())
          return {8.49e-01, 1.15e+00};
        else if (type == residual_estimator_ESV07_id())
          return {0.0, 0.0};
        else if (type == diffusive_flux_estimator_id())
          return {4.76e-02, 3.06e-02};
        else if (type == efficiency_ESV07_id())
          return {2.57e+01, 6.99e+01};
        else
          return BaseType::expected_results(type);
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
#endif
      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  } // ... expected_results(...)

  virtual std::map< std::string, std::vector< double > > run(std::ostream& out = std::cout)
  {
    return BaseType::BaseType::run(false, out);
  }

private:
  double compute_energy_norm()
  {
    using namespace Dune;
    using namespace GDT;
    // get current solution
    assert(current_level_ < test_.num_levels());
    if (last_computed_level_ != current_level_) {
      this->compute_on_current_refinement();
    }
    assert(last_computed_level_ == current_level_);
    const DiscretizationType discretization(test_.level_grid_part(current_level_), test_.boundary_info(),
                                            test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
    assert(current_solution_vector_on_level_);
    const ConstDiscreteFunctionType current_solution(discretization.space(),
                                                      *current_solution_vector_on_level_,
                                                      "discrete solution");
    // compute error
    if (test_.provides_exact_solution()) {
      typedef Dune::Stuff::Functions::Difference< ExactSolutionType, ConstDiscreteFunctionType > DifferenceType;
      const DifferenceType difference(test_.exact_solution(), current_solution);
      auto grid_view = test_.level_grid_view(current_level_);
      const GDT::Products::Elliptic< GridViewType, typename TestCase::DiffusionType >
          elliptic_product(grid_view, test_.diffusion(), over_integrate);
      return std::sqrt(elliptic_product.apply2(difference, difference));
    } else {
      if (!reference_solution_computed_)
        this->compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      assert(current_solution_vector_);
      const VectorType difference_vector = (*reference_solution_vector_) - (*current_solution_vector_);
      const ConstDiscreteFunctionType difference(reference_discretization_->space(), difference_vector);
      auto grid_view = test_.reference_grid_view();
      const GDT::Products::Elliptic< GridViewType, typename TestCase::DiffusionType >
          elliptic_product(grid_view, test_.diffusion(), over_integrate);
      return std::sqrt(elliptic_product.apply2(difference, difference));
    }
  } // ... compute_energy_norm(...)

  double compute_nonconformity_estimator()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    BaseType::compute_on_current_refinement();
    assert(current_solution_vector_on_level_);
    const auto grid_part = test_.level_grid_part(current_level_);
    auto grid_view = test_.level_grid_view(current_level_);
    const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                            test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(discretization.space(), *current_solution_vector_on_level_);
    VectorType oswald_interpolation_vector(discretization.space().mapper().size());
    DiscreteFunctionType oswald_interpolation(discretization.space(), oswald_interpolation_vector);

    const Operators::OswaldInterpolation< GridViewType > oswald_interpolation_operator(grid_view);
    oswald_interpolation_operator.apply(discrete_solution, oswald_interpolation);

    const Products::Elliptic< GridViewType, typename TestCase::DiffusionType >
        elliptic_product(grid_view, test_.diffusion(), over_integrate);
    return elliptic_product.induced_norm(discrete_solution - oswald_interpolation);
  } // ... compute_nonconformity_estimator(...)

  double compute_residual_estimator_ESV07()
  {
    using namespace Dune;
    using namespace GDT;
    const auto grid_view = test_.level_grid_view(current_level_);

    typedef Spaces::FV::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
    const P0SpaceType p0_space(grid_view);
    typedef DiscreteFunction< P0SpaceType, VectorType > P0DiscreteFunctionType;
    P0DiscreteFunctionType p0_force(p0_space);

    Operators::Projection< GridViewType > projection_operator(grid_view, over_integrate);
    projection_operator.apply(test_.force(), p0_force);

    typedef typename Stuff::Functions::ESV2007::Cutoff< typename TestCase::DiffusionType > CutoffFunctionType;
    const CutoffFunctionType cutoff_function(test_.diffusion());

    const Products::WeightedL2< GridViewType, CutoffFunctionType >
        weighted_l2_product(grid_view, cutoff_function, over_integrate);
    return weighted_l2_product.induced_norm(test_.force() - p0_force);
  } // ... compute_residual_estimator_ESV07(...)

  double compute_diffusive_flux_estimator()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    BaseType::compute_on_current_refinement();
    assert(current_solution_vector_on_level_);

    const auto grid_part = test_.level_grid_part(current_level_);
    const auto grid_view = test_.level_grid_view(current_level_);
    const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                            test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(discretization.space(), *current_solution_vector_on_level_);

    typedef Spaces::RT::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
    const RTN0SpaceType rtn0_space(grid_view);
    typedef DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
    RTN0DiscreteFunctionType diffusive_flux(rtn0_space);

    typedef typename TestCase::DiffusionType DiffusionType;
    const Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionType >
      diffusive_flux_reconstruction(grid_view, test_.diffusion(), over_integrate);
    diffusive_flux_reconstruction.apply(discrete_solution, diffusive_flux);

    GDT::Products::ESV2007::DiffusiveFluxEstimate< GridViewType, DiffusionType, RTN0DiscreteFunctionType
                                                , ConstDiscreteFunctionType, ConstDiscreteFunctionType >
      diffusive_flux_estimator_product(grid_view, discrete_solution, discrete_solution,
                                       test_.diffusion(), diffusive_flux, over_integrate);
    return std::sqrt(diffusive_flux_estimator_product.apply2());
  } // ... compute_diffusive_flux_estimator(...)

  double compute_estimator_ESV07()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    BaseType::compute_on_current_refinement();
    assert(current_solution_vector_on_level_);
    const auto grid_part = test_.level_grid_part(current_level_);
    const auto grid_view = test_.level_grid_view(current_level_);

    const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                            test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(discretization.space(), *current_solution_vector_on_level_);

    VectorType oswald_interpolation_vector(discretization.space().mapper().size());
    DiscreteFunctionType oswald_interpolation(discretization.space(), oswald_interpolation_vector);
    const Operators::OswaldInterpolation< GridViewType > oswald_interpolation_operator(grid_view);
    oswald_interpolation_operator.apply(discrete_solution, oswald_interpolation);

    typedef Spaces::FV::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
    const P0SpaceType p0_space(grid_view);
    VectorType p0_force_vector(p0_space.mapper().size());
    typedef DiscreteFunction< P0SpaceType, VectorType > P0DiscreteFunctionType;
    P0DiscreteFunctionType p0_force(p0_space, p0_force_vector);
    Operators::Projection< GridViewType > projection_operator(grid_view);
    projection_operator.apply(test_.force(), p0_force);

    typedef Spaces::RT::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
    const RTN0SpaceType rtn0_space(grid_view);
    VectorType diffusive_flux_vector(rtn0_space.mapper().size());
    typedef DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
    RTN0DiscreteFunctionType diffusive_flux(rtn0_space, diffusive_flux_vector);

    typedef typename TestCase::DiffusionType DiffusionType;
    const Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionType >
      diffusive_flux_reconstruction(grid_view, test_.diffusion());
    diffusive_flux_reconstruction.apply(discrete_solution, diffusive_flux);

    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionType > >
      local_eta_nc_product(1, test_.diffusion());
    const auto eta_nc_difference = discrete_solution - oswald_interpolation;

    typedef typename Stuff::Functions::ESV2007::Cutoff< DiffusionType > CutoffFunctionType;
    const CutoffFunctionType cutoff_function(test_.diffusion());
    const  LocalOperator::Codim0Integral< LocalEvaluation::Product< CutoffFunctionType > >
      local_eta_r_product(1, cutoff_function);
    const auto eta_r_difference = test_.force() - p0_force;

    const LocalOperator::Codim0Integral< LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionType
                                                                                        , RTN0DiscreteFunctionType > >
      local_eta_df_product(1, test_.diffusion(), diffusive_flux);

    // walk the grid
    double eta = 0.0;
    std::vector< DynamicMatrix< RangeFieldType > > tmp_matrices(std::max(std::max(local_eta_nc_product.numTmpObjectsRequired(),
                                                                                  local_eta_r_product.numTmpObjectsRequired()),
                                                                         local_eta_df_product.numTmpObjectsRequired()),
                                                                DynamicMatrix< RangeFieldType >(1, 1, 0.0));
    DynamicMatrix< RangeFieldType > local_result_matrix(1, 1, 0.0);
    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;

      local_result_matrix *= 0.0;
      const auto local_eta_nc_difference = eta_nc_difference.local_function(entity);
      local_eta_nc_product.apply(*local_eta_nc_difference, *local_eta_nc_difference, local_result_matrix, tmp_matrices);
      assert(local_result_matrix.rows() >= 1);
      assert(local_result_matrix.cols() >= 1);
      const double eta_nc_t_squared = local_result_matrix[0][0];

      local_result_matrix *= 0.0;
      const auto local_eta_r_difference = eta_r_difference.local_function(entity);
      local_eta_r_product.apply(*local_eta_r_difference, *local_eta_r_difference, local_result_matrix, tmp_matrices);
      assert(local_result_matrix.rows() >= 1);
      assert(local_result_matrix.cols() >= 1);
      const double eta_r = std::sqrt(local_result_matrix[0][0]);

      local_result_matrix *= 0.0;
      const auto local_discrete_solution = discrete_solution.local_function(entity);
      local_eta_df_product.apply(*local_discrete_solution, *local_discrete_solution, local_result_matrix, tmp_matrices);
      assert(local_result_matrix.rows() >= 1);
      assert(local_result_matrix.cols() >= 1);
      const double eta_df = std::sqrt(local_result_matrix[0][0]);

      eta += eta_nc_t_squared + std::pow(eta_r + eta_df, 2);
    } // walk the grid
    return std::sqrt(eta);
  } // ... compute_estimator_ESV07(...)

  double compute_efficiency_ESV07()
  {
    return compute_estimator_ESV07() / compute_energy_norm();
  }

#if 0
  double compute_residual_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const RangeFieldType poincare_constant_sqrt = RangeFieldType(1) / M_PIl;
    const size_t integration_order = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                            test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(discretization.space(),
                                                      *current_solution_vector_on_level_,
                                                      "discrete solution");

    typedef RaviartThomasSpace::FemLocalfunctionsWrapper< GridPartType, 0, RangeFieldType, dimDomain > RTN_SpaceType;
    const RTN_SpaceType rtn_space(grid_part);
    typedef DiscreteFunction< RTN_SpaceType, VectorType > RTN_DiscreteFunctionType;
    VectorType rtn_vector(rtn_space.mapper().size());
    RTN_DiscreteFunctionType diffusive_flux_reconstruction(rtn_space,
                                                           rtn_vector,
                                                           "diffusive flux reconstruction");
    std::vector< typename RTN_SpaceType::BaseFunctionSetType::RangeType >
        rtn_basis_values(rtn_space.mapper().maxNumDofs(),
                         typename RTN_SpaceType::BaseFunctionSetType::RangeType(0));

    // walk the grid for the second time
    const auto entity_it_end = grid_part->template end< 0 >();
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_solution_entity = discrete_solution.local_function(entity);
      const auto local_diffusion_entity = test_.diffusion().local_function(entity);
      auto local_diffusive_flux = diffusive_flux_reconstruction.local_discrete_function(entity);
      auto local_diffusive_flux_DoF_vector = local_diffusive_flux.vector();
      const auto rtn_basis = rtn_space.baseFunctionSet(entity);
      // get the local finite elements
      const auto rtn_finite_element = rtn_space.backend().finiteElement(entity);
      const auto& rtn_local_coefficients = rtn_finite_element.localCoefficients();
      std::vector< size_t > intersection_id_to_local_DoF_id_map(rtn_local_coefficients.size(), 0);
      for (size_t ii = 0; ii < rtn_local_coefficients.size(); ++ii)
        intersection_id_to_local_DoF_id_map[rtn_local_coefficients.localKey(ii).subEntity()] = ii;
      // to compute the diffusive flux reconstruction
      // * loop over all intersections
      const auto intersection_end_it = grid_part->iend(entity);
      for (auto intersection_it = grid_part->ibegin(entity);
           intersection_it != intersection_end_it;
           ++intersection_it) {
        // get intersection information
        const auto& intersection = *intersection_it;
        const size_t intersection_index = intersection.indexInInside();
        assert(intersection_index < intersection_id_to_local_DoF_id_map.size() && "This should not happen!");
        const size_t intersection_DoF_index = intersection_id_to_local_DoF_id_map[intersection_index];
        // prepare quadrature
        const size_t quadrature_order = 1;
        const auto& face_quadrature = QuadratureRules< DomainFieldType, dimDomain - 1 >::rule(intersection.type(),
                                                                                              2*quadrature_order + 1);
        RangeFieldType lhs_integral = 0;
        RangeFieldType rhs_integral = 0;

        if (intersection.boundary() && !intersection.neighbor()) {
          // we are on the domain boundary
          // loop over all quadrature points
          for (auto quadrature_point : face_quadrature) {
            const FieldVector< DomainFieldType, dimDomain - 1 >& point_intersection = quadrature_point.position();
            const FieldVector< DomainFieldType, dimDomain >
                point_entity = intersection.geometryInInside().global(point_intersection);
            const RangeFieldType integration_factor = intersection.geometry().integrationElement(point_intersection);
            const RangeFieldType quadrature_weight = quadrature_point.weight();
            // evaluate
            const auto unit_outer_normal = intersection.unitOuterNormal(point_intersection);
            rtn_basis.evaluate(point_entity, rtn_basis_values);
            const auto& rtn_basis_value = rtn_basis_values[intersection_DoF_index];
            const auto solution_value = local_solution_entity->evaluate(point_entity);
            const auto solution_gradient = local_solution_entity->jacobian(point_entity);
            const auto diffusion_value = local_diffusion_entity->evaluate(point_entity);
            const RangeFieldType sigma = 14.0;
            const RangeFieldType gamma = diffusion_value;
            const RangeFieldType penalty = (sigma * gamma ) / std::pow(intersection.geometry().volume(), 1.0);
            // compute integrals
            lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
            rhs_integral += integration_factor * quadrature_weight
                            * (-1.0 * diffusion_value * (solution_gradient[0] * unit_outer_normal)
                                + penalty * solution_value);
          } // loop over all quadrature points
          // set local DoF
          local_diffusive_flux_DoF_vector.set(intersection_DoF_index, rhs_integral / lhs_integral);
        } else if (intersection.neighbor() && !intersection.boundary()) {
          // we are on the inside
          // get the neighbour
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // work only once on each intersection
          if (grid_part->indexSet().index(entity) < grid_part->indexSet().index(neighbour)) {
            // get the neighbours local functions
            const auto local_solution_neighbour = discrete_solution.local_function(neighbour);
            const auto local_diffusion_neighbour = test_.diffusion().local_function(neighbour);
            // loop over all quadrature points
            for (auto quadrature_point : face_quadrature) {
              const FieldVector< DomainFieldType, dimDomain - 1 > point_intersection = quadrature_point.position();
              const FieldVector< DomainFieldType, dimDomain >
                  point_entity = intersection.geometryInInside().global(point_intersection);
              const FieldVector< DomainFieldType, dimDomain >
                  point_neighbour = intersection.geometryInOutside().global(point_intersection);
              const RangeFieldType integration_factor = intersection.geometry().integrationElement(point_intersection);
              const RangeFieldType quadrature_weight = quadrature_point.weight();
              // evaluate
              const auto unit_outer_normal = intersection.unitOuterNormal(point_intersection);
              rtn_basis.evaluate(point_entity, rtn_basis_values);
              const auto& rtn_basis_value = rtn_basis_values[intersection_DoF_index];
              const auto solution_value_entity = local_solution_entity->evaluate(point_entity);
              const auto solution_value_neighbour = local_solution_entity->evaluate(point_neighbour);
              const auto solution_gradient_entity = local_solution_entity->jacobian(point_entity);
              const auto solution_gradient_neighbour = local_solution_neighbour->jacobian(point_neighbour);
              const auto diffusion_value_entity = local_diffusion_entity->evaluate(point_entity);
              const auto diffusion_value_neighbour = local_diffusion_neighbour->evaluate(point_neighbour);
              // evaluate penalty parameter and weighting
              const RangeFieldType sigma = 8.0;
              const RangeFieldType delta_plus  = diffusion_value_neighbour;
              const RangeFieldType delta_minus = diffusion_value_entity;
              const RangeFieldType gamma = (delta_plus * delta_minus)/(delta_plus + delta_minus);
              const RangeFieldType penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), 1.0);
              const RangeFieldType weight_plus = delta_minus / (delta_plus + delta_minus);
              const RangeFieldType weight_minus = delta_plus / (delta_plus + delta_minus);
              // compute integrals
              lhs_integral += integration_factor * quadrature_weight * (rtn_basis_value * unit_outer_normal);
              rhs_integral += integration_factor * quadrature_weight
                              * ( - weight_minus * (diffusion_value_entity
                                                        * (solution_gradient_entity[0] * unit_outer_normal))
                                  - weight_plus * (diffusion_value_neighbour
                                                          * (solution_gradient_neighbour[0] * unit_outer_normal))
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
    DiscreteFunction< RTN_SpaceType, VectorType > diffusive_flux_reconstruction_2(rtn_space,
                                                                                diffusive_flux_reconstruction_vector_2,
                                                                                "diffusive flux reconstruction");
    // reconstruct
    typedef Operators::DiffusiveFluxReconstruction< GridPartType,
                                                   typename TestCase::DiffusionType > ReconstructionOperatorType;
    const ReconstructionOperatorType reconstruction_operator(*grid_part, test_.diffusion());
    reconstruction_operator.apply(discrete_solution, diffusive_flux_reconstruction_2);
//    const auto difference = diffusive_flux_reconstruction.vector().backend() - diffusive_flux_reconstruction_2.vector().backend();
//    const double max = std::max(std::abs(difference.maxCoeff()), std::abs(difference.minCoeff()));
//    if (max > 0.0)
//      DUNE_THROW(InvalidStateException, "boom " << max);

    // walk the grid for the third time
    std::vector< RangeFieldType > estimators_residual(grid_part->indexSet().size(0), RangeFieldType(0));
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion = test_.diffusion().local_function(entity);
      const auto local_force = test_.force().local_function(entity);
      const auto local_diffusive_flux_reconstruction = diffusive_flux_reconstruction_2.local_function(entity);
      RangeFieldType diffusion_min = std::numeric_limits< RangeFieldType >::max();
      // do a volume quadrature
      const auto& volume_quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                          2*integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector< DomainFieldType, dimDomain >& point_entity = quadrature_point.position();
        const double quadrature_weight = quadrature_point.weight();
        // evaluate
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value = local_diffusion->evaluate(point_entity);
        if (diffusion_value < diffusion_min)
          diffusion_min = diffusion_value;
        const auto force_value = local_force->evaluate(point_entity);
        const auto diffusive_flux_gradient = local_diffusive_flux_reconstruction->jacobian(point_entity);
        // compute the residual estimator
        const auto diffusive_flux_divergence = diffusive_flux_gradient[0][0] + diffusive_flux_gradient[1][1];
        const auto residual_difference = force_value - diffusive_flux_divergence;
        const auto residual_product = residual_difference * residual_difference;
        estimators_residual[entity_index] += integration_factor * quadrature_weight * residual_product;
      } // do a volume quadrature
      // for the residual estimator
      estimators_residual[entity_index] *= (entity.geometry().volume() * poincare_constant_sqrt)
                                           / std::sqrt(diffusion_min);
    } // walk the grid for the third time

    // compute error components
    double residual_estimator = 0;
    for (size_t ii = 0; ii < estimators_residual.size(); ++ii) {
      residual_estimator += std::pow(estimators_residual[ii], 2);
    }
    return std::sqrt(residual_estimator);
  } // ... compute_residual_estimator_ESV10(...)
#endif

#if 0
  double compute_estimator_ESV10()
  {
    using namespace Dune;
    using namespace Dune::GDT;

    const auto poincare_constant_sqrt = RangeFieldType(1) / M_PIl;
    const size_t integration_order = 5;

    // prepare discrete solution
    BaseType::compute_on_current_refinement();
    const auto grid_part = test_.level_grid_part(current_level_);
    const DiscretizationType discretization(grid_part, test_.boundary_info(), test_.diffusion(), test_.force(),
                                            test_.dirichlet(), test_.neumann());
    const ConstDiscreteFunctionType discrete_solution(discretization.space(),
                                                      *current_solution_vector_on_level_,
                                                      "discrete solution");
    VectorType oswald_projection_vector(discretization.space().mapper().size());
    DiscreteFunctionType oswald_projection(discretization.space(),
                                           oswald_projection_vector,
                                           "oswald projection");

    typedef typename DiscretizationType::SpaceType TestSpaceType;
    typedef FieldVector< DomainFieldType, dimDomain > DomainType;

    typedef RaviartThomasSpace::FemLocalfunctionsWrapper< GridPartType, 0, RangeFieldType, dimDomain > RTN_SpaceType;
    const RTN_SpaceType rtn_space(grid_part);
    typedef DiscreteFunction< RTN_SpaceType, VectorType > RTN_DiscreteFunctionType;
    VectorType rtn_vector(rtn_space.mapper().size());
    RTN_DiscreteFunctionType diffusive_flux_reconstruction(rtn_space,
                                                           rtn_vector,
                                                           "diffusive flux reconstruction");
    // reconstruct
    typedef Operators::DiffusiveFluxReconstruction< GridPartType,
                                                   typename TestCase::DiffusionType > ReconstructionOperatorType;
    const ReconstructionOperatorType reconstruction_operator(*grid_part, test_.diffusion());
    reconstruction_operator.apply(discrete_solution, diffusive_flux_reconstruction);

    // data structures we need
    // * a map from a global vertex id to global DoF ids
    //   given a vertex, one obtains a set of all global DoF ids, which are associated with this vertex
    typedef std::vector< std::set< size_t > > VertexToEntitiesMapType;
    VertexToEntitiesMapType vertex_to_dof_id_map(grid_part->indexSet().size(dimDomain));
    // * a set to hold the global id off all boundary vertices
    std::set< size_t > boundary_vertices;
    // * vectors to hold the local estimators
    std::vector< RangeFieldType > estimators_nonconformity(grid_part->indexSet().size(0), RangeFieldType(0));
    std::vector< RangeFieldType > estimators_residual(grid_part->indexSet().size(0), RangeFieldType(0));
    std::vector< RangeFieldType > estimators_diffusive_flux(grid_part->indexSet().size(0), RangeFieldType(0));

    // walk the grid for the first time
    const auto entity_it_end = grid_part->template end< 0 >();
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local finite elements
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size()
             && "Wrong finite element given!");
      // loop over all vertices
      std::vector< DomainType > global_vertices(entity.template count< dimDomain >(), DomainType(0));
      std::vector< size_t > global_vertex_ids(global_vertices.size(), 0);
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get global vertex id
        const auto vertexPtr = entity.template subEntity< dimDomain >(local_vertex_id);
        const auto& vertex = *vertexPtr;
        global_vertex_ids[local_vertex_id] = grid_part->indexSet().index(vertex);
        global_vertices[local_vertex_id] = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            const size_t local_DOF_id = entity_dg_local_key.index();
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
               local_intersection_corner_id < boost::numeric_cast< size_t >(intersection_geometry.corners());
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
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_solution_entity = discrete_solution.local_discrete_function(entity);
      const auto& local_solution_entity_DoF_vector = local_solution_entity.vector();
      // get the local finite elements
      // * for the oswald projection
      typedef typename TestSpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = discretization.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size()
             && "Wrong finite element given!");
      // * for the diffusive flux reconstruction
      const auto rtn_finite_element = rtn_space.backend().finiteElement(entity);
      const auto& rtn_local_coefficients = rtn_finite_element.localCoefficients();
      std::vector< size_t > intersection_id_to_local_DoF_id_map(rtn_local_coefficients.size(), 0);
      for (size_t ii = 0; ii < rtn_local_coefficients.size(); ++ii)
        intersection_id_to_local_DoF_id_map[rtn_local_coefficients.localKey(ii).subEntity()] = ii;
      // to compute the oswald projection
      // * loop over all local DoFs
      for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
        const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
        const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
        assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
        const size_t local_vertex_id = entity_cg_local_key.subEntity();
        const size_t local_DoF_id = entity_dg_local_key.index();
        const auto vertexPtr = entity.template subEntity< dimDomain >(local_vertex_id);
        const auto& vertex = *vertexPtr;
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
    for (auto entity_it = grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const size_t entity_index = grid_part->indexSet().index(entity);
      // get the local functions
      const auto local_diffusion = test_.diffusion().local_function(entity);
      const auto local_force = test_.force().local_function(entity);
      const auto local_solution = discrete_solution.local_function(entity);
      const auto local_oswald_projection = oswald_projection.local_function(entity);
      const auto local_diffusive_flux_reconstruction = diffusive_flux_reconstruction.local_function(entity);
      RangeFieldType diffusion_min = std::numeric_limits< RangeFieldType >::max();
      // do a volume quadrature
      const auto& volume_quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(entity.type(),
                                                                                          2*integration_order + 1);
      for (auto quadrature_point : volume_quadrature) {
        const FieldVector< DomainFieldType, dimDomain > point_entity = quadrature_point.position();
        const double quadrature_weight = quadrature_point.weight();
        // evaluate
        const double integration_factor = entity.geometry().integrationElement(point_entity);
        const auto diffusion_value = local_diffusion->evaluate(point_entity);
        if (diffusion_value < diffusion_min)
          diffusion_min = diffusion_value;
        const auto force_value = local_force->evaluate(point_entity);
        auto solution_gradient = local_solution->jacobian(point_entity);
        const auto oswald_projection_gradient = local_oswald_projection->jacobian(point_entity);
        auto diffusive_flux_value = local_diffusive_flux_reconstruction->evaluate(point_entity);
        const auto diffusive_flux_gradient = local_diffusive_flux_reconstruction->jacobian(point_entity);
        // compute local nonconformity estimator
        const auto nonconformity_difference = solution_gradient[0] - oswald_projection_gradient[0];
        const auto nonconformity_product = nonconformity_difference * nonconformity_difference;
        estimators_nonconformity[entity_index] += integration_factor * quadrature_weight
                                                  * diffusion_value * nonconformity_product;
        // compute the residual estimator
        const auto diffusive_flux_divergence = diffusive_flux_gradient[0][0] + diffusive_flux_gradient[1][1];
        const auto residual_difference = force_value - diffusive_flux_divergence;
        const auto residual_product = residual_difference * residual_difference;
        estimators_residual[entity_index] += integration_factor * quadrature_weight * residual_product;
        // compute diffusive flux estimator
        const auto diffusion_value_sqrt = std::sqrt(diffusion_value);
        solution_gradient[0] *= diffusion_value_sqrt;
        diffusive_flux_value /= diffusion_value_sqrt;
        const auto diffusive_flux_sum = solution_gradient[0] + diffusive_flux_value;
        const auto diffusive_flux_product = diffusive_flux_sum * diffusive_flux_sum;
        estimators_diffusive_flux[entity_index] += integration_factor * quadrature_weight * diffusive_flux_product;
      } // do a volume quadrature
      // for the residual estimator
      estimators_residual[entity_index] *= (entity.geometry().volume() * poincare_constant_sqrt)
                                           / std::sqrt(diffusion_min);
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
#endif

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
