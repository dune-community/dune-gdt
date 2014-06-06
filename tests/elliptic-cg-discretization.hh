// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH
#define DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH

#if !HAVE_DUNE_PDELAB
#error "This one requires dune-pdelab!"
#endif

#include <memory>
#include <vector>
#include <limits>
#include <type_traits>

#include <dune/common/timer.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/convergence-study.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/spaces/continuouslagrange/pdelab.hh>
#include <dune/gdt/operators/elliptic-cg.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/spaces/constraints.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/operators/prolongations.hh>
#include <dune/gdt/operators/projections.hh>

#include "elliptic-testcases.hh"


namespace EllipticCG {


template <class GridViewType, int polynomialOrder,
          class MatrixImp = typename Dune::Stuff::LA::Container<double>::MatrixType,
          class VectorImp = typename Dune::Stuff::LA::Container<double>::VectorType>
class Discretization
{
public:
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef typename GridViewType::ctype DomainFieldType;

  static const unsigned int dimRange = 1;
  typedef double RangeFieldType;

  static const unsigned int polOrder = polynomialOrder;

  typedef Dune::Stuff::Grid::BoundaryInfoInterface<typename GridViewType::Intersection> BoundaryInfoType;
  typedef Dune::Stuff::LocalizableFunctionInterface<typename GridViewType::template Codim<0>::Entity, DomainFieldType,
                                                    dimDomain, RangeFieldType, dimRange> FunctionType;

  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;

  typedef Dune::GDT::Spaces::ContinuousLagrange::PdelabBased<GridViewType, polOrder, RangeFieldType, dimRange>
      SpaceType;

  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::ConstDiscreteFunction<SpaceType, VectorType> ConstDiscreteFunctionType;

  Discretization(const std::shared_ptr<const GridViewType>& gp, const BoundaryInfoType& info, const FunctionType& diff,
                 const FunctionType& forc, const FunctionType& dir, const FunctionType& neu)
    : space_(gp)
    , boundary_info_(info)
    , diffusion_(diff)
    , force_(forc)
    , dirichlet_(dir)
    , neumann_(neu)
    , is_assembled_(false)
    , system_matrix_(0, 0)
    , rhs_vector_(0)
    , dirichlet_shift_vector_(0)
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
    using namespace Dune;
    using namespace Dune::GDT;
    if (!is_assembled_) {
      // create the containers (use the sparsity pattern of the operator)
      typedef Operators::EllipticCG<FunctionType, MatrixType, SpaceType> EllipticOperatorType;
      system_matrix_ =
          MatrixType(space_.mapper().size(), space_.mapper().size(), EllipticOperatorType::pattern(space_));
      rhs_vector_             = VectorType(space_.mapper().size());
      dirichlet_shift_vector_ = VectorType(space_.mapper().size());

      // define the lhs operator and the rhs functionals
      EllipticOperatorType elliptic_operator(diffusion_, system_matrix_, space_);
      typedef GDT::Functionals::L2Volume<FunctionType, VectorType, SpaceType> L2VolumeFunctionalType;
      L2VolumeFunctionalType force_functional(force_, rhs_vector_, space_);
      typedef GDT::Functionals::L2Face<FunctionType, VectorType, SpaceType> L2FaceFunctionalType;
      L2FaceFunctionalType neumann_functional(neumann_, rhs_vector_, space_);
      // project the dirichlet boundary values
      DiscreteFunctionType dirichlet_projection(space_, dirichlet_shift_vector_);
      typedef Operators::DirichletProjectionLocalizable<GridViewType, FunctionType, DiscreteFunctionType>
          DirichletProjectionOperator;
      DirichletProjectionOperator dirichlet_projection_operator(
          *(space_.grid_view()), boundary_info_, dirichlet_, dirichlet_projection);
      // assemble everything
      SystemAssembler<SpaceType> grid_walker(space_);
      grid_walker.add(elliptic_operator);
      grid_walker.add(force_functional);
      grid_walker.add(neumann_functional, new ApplyOn::NeumannIntersections<GridViewType>(boundary_info_));
      grid_walker.add(dirichlet_projection_operator, new ApplyOn::BoundaryEntities<GridViewType>());
      grid_walker.walk();
      // substract the operators action on the dirichlet values
      auto tmp = rhs_vector_.copy();
      system_matrix_.mv(dirichlet_shift_vector_, tmp);
      rhs_vector_ -= tmp;
      // apply the dirichlet constraints
      Constraints::Dirichlet<typename GridViewType::Intersection, RangeFieldType> dirichlet_constraints(
          boundary_info_, space_.mapper().maxNumDofs(), space_.mapper().maxNumDofs());
      grid_walker.add(dirichlet_constraints, system_matrix_);
      grid_walker.add(dirichlet_constraints, rhs_vector_);
      grid_walker.walk();

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

  const VectorType& dirichlet_shift_vector() const
  {
    return dirichlet_shift_vector_;
  }

  void solve(VectorType& solution) const
  {
    if (!is_assembled_)
      assemble();

    // solve
    Dune::Stuff::LA::Solver<MatrixType>(system_matrix_).apply(rhs_vector_, solution);
    solution += dirichlet_shift_vector_;
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
  mutable bool is_assembled_;
  mutable MatrixType system_matrix_;
  mutable VectorType rhs_vector_;
  mutable VectorType dirichlet_shift_vector_;
}; // class Discretization


template <class TestCase, int polOrder>
class EocStudy : public Dune::Stuff::Common::ConvergenceStudy
{
  typedef Dune::Stuff::Common::ConvergenceStudy BaseType;

protected:
  typedef typename TestCase::GridViewType GridViewType;
  typedef typename TestCase::EntityType EntityType;

  typedef typename TestCase::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = TestCase::dimDomain;
  typedef typename TestCase::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = TestCase::dimRange;

  typedef Discretization<GridViewType, polOrder> DiscretizationType;
  typedef typename DiscretizationType::VectorType VectorType;
  typedef typename DiscretizationType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscretizationType::ConstDiscreteFunctionType ConstDiscreteFunctionType;

  typedef typename TestCase::ExactSolutionType ExactSolutionType;

public:
  EocStudy(const TestCase& test)
    : test_(test)
    , current_level_(0)
    , last_computed_level_(std::numeric_limits<size_t>::max())
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
    return "continuous galerkin discretization, polOrder " + Dune::Stuff::Common::toString(polOrder);
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
          reference_discretization_->space(), *reference_solution_vector_, "CG reference solution");
      // compute norm
      return compute_norm(*(test_.reference_grid_view()), reference_solution, type);
    }
  } // ... norm_reference_solution(...)

  virtual size_t current_grid_size() const DS_OVERRIDE DS_FINAL
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_view = test_.level_grid_view(current_level_);
    return grid_view->size(0);
  } // ... current_grid_size(...)

  virtual double current_grid_width() const DS_OVERRIDE DS_FINAL
  {
    assert(current_level_ < test_.num_levels());
    const auto grid_view  = test_.level_grid_view(current_level_);
    double h              = 0.0;
    const auto entity_ptr = grid_view->template begin<0>();
    const auto& entity = *entity_ptr;
    for (int cc = 0; cc < entity.template count<dimDomain>(); ++cc) {
      const auto vertex = entity.template subEntity<dimDomain>(cc)->geometry().center();
      for (int dd = cc + 1; dd < entity.template count<dimDomain>(); ++dd) {
        const auto other_vertex = entity.template subEntity<dimDomain>(dd)->geometry().center();
        const auto diff         = vertex - other_vertex;
        h                       = std::max(h, diff.two_norm());
      }
    }
    return h;
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
      const auto grid_view = test_.level_grid_view(current_level_);
      const DiscretizationType discretization(
          grid_view, test_.boundary_info(), test_.diffusion(), test_.force(), test_.dirichlet(), test_.neumann());
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
      const Operators::Prolongation<GridViewType> prolongation_operator(*reference_grid_view);
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
      return compute_norm(*(test_.reference_grid_view()), test_.exact_solution() - current_solution, type);
    } else {
      // get reference solution
      compute_reference_solution();
      assert(reference_discretization_);
      assert(reference_solution_vector_);
      const ConstDiscreteFunctionType reference_solution(
          reference_discretization_->space(), *reference_solution_vector_, "CG reference solution");
      return compute_norm(*(test_.reference_grid_view()), reference_solution - current_solution, type);
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
          return {1.97e-01, 4.86e-02, 1.22e-02, 3.03e-03};
        else if (type.compare("H1_semi") == 0)
          return {4.12e-01, 2.08e-01, 1.04e-01, 5.18e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::LocalThermalBlock<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {7.93e-02, 4.16e-02, 1.20e-02, 2.73e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.52e-01, 3.03e-01, 1.64e-01, 7.52e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::ER07<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {1.49e-01, 3.83e-02, 9.66e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.60e-01, 1.85e-01, 9.25e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::MixedBoundaryTypes<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {8.32e-02, 2.23e-02, 5.53e-03, 1.20e-03};
        else if (type.compare("H1_semi") == 0)
          return {3.12e-01, 1.65e-01, 8.24e-02, 3.76e-02};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else if (std::is_same<TestCase, EllipticTestCase::Spe10Model1<Dune::ALUConformGrid<2, 2>>>::value) {
      if (polOrder == 1) {
        if (type.compare("L2") == 0)
          return {7.80e-02, 4.91e-02};
        else if (type.compare("H1_semi") == 0)
          return {4.04e-01, 3.78e-01};
        else
          DUNE_THROW(Dune::RangeError, "Wrong type '" << type << "' requested!");
      } else
        DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this polOrder!");
    } else
      DUNE_THROW(Dune::NotImplemented, "Please record the expected results for this TestCase/GridType combination!");
  } // ... expected_results(...)

  virtual std::map<std::string, std::vector<double>> run(std::ostream& out = std::cout)
  {
    return BaseType::run(true, out);
  }

private:
  void compute_reference_solution()
  {
    if (!reference_solution_computed_) {
      reference_discretization_ =
          std::unique_ptr<DiscretizationType>(new DiscretizationType(test_.reference_grid_view(),
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

  template <class GridViewType, class FunctionType>
  double compute_norm(const GridViewType& grid_view, const FunctionType& function, const std::string type)
  {
    using namespace Dune;
    using namespace Dune::GDT;
    if (type.compare("L2") == 0) {
      Products::L2<GridViewType> l2_product_operator(grid_view);
      return std::sqrt(l2_product_operator.apply2(function, function));
    } else if (type.compare("H1_semi") == 0) {
      Products::H1SemiGeneric<GridViewType> h1_product_operator(grid_view);
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


} // namespace EllipticCG

#endif // DUNE_GDT_TEST_ELLIPTIC_CG_DISCRETIZATION_HH
