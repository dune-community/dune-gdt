// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH

#include <utility>
#include <vector>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/timestepper/explicit-rungekutta.hh>
#include <dune/gdt/timestepper/fractional-step.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Discretizations {


// forward
template< class ProblemType, class AnsatzSpaceType, class MatrixType, class VectorType, class TestSpaceType = AnsatzSpaceType >
class StationaryContainerBasedDefault;

template< class ProblemImp, class FVSpaceImp >
class NonStationaryDefault;


namespace internal {


template< class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp >
class StationaryContainerBasedDefaultTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef StationaryContainerBasedDefault
      < ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp > derived_type;
  typedef ProblemImp     ProblemType;
  typedef AnsatzSpaceImp AnsatzSpaceType;
  typedef TestSpaceImp   TestSpaceType;
  typedef MatrixImp      MatrixType;
  typedef VectorImp      VectorType;
}; // class StationaryContainerBasedDefaultTraits


template< class ProblemImp, class FVSpaceImp >
class NonStationaryDefaultTraits
{
  // no checks of the arguments needed, those are done in the interfaces
public:
  typedef NonStationaryDefault
      < ProblemImp, FVSpaceImp >                                          derived_type;
  typedef ProblemImp                                                      ProblemType;
  typedef FVSpaceImp                                                      FVSpaceType;
  typedef typename FVSpaceType::RangeFieldType                            RangeFieldType;
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType >   VectorType;
  typedef DiscreteFunction< FVSpaceType, VectorType >                     DiscreteFunctionType;
  typedef std::map< double, DiscreteFunctionType, Dune::GDT::TimeStepper::internal::FloatCmpLt > DiscreteSolutionType;
}; // class NonStationaryDefaultTraits


} // namespace internal


template< class ProblemImp, class AnsatzSpaceImp, class MatrixImp, class VectorImp, class TestSpaceImp >
class StationaryContainerBasedDefault
  : public ContainerBasedStationaryDiscretizationInterface<
             internal::StationaryContainerBasedDefaultTraits< ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp > >
{
  typedef ContainerBasedStationaryDiscretizationInterface
          < internal::StationaryContainerBasedDefaultTraits< ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp > >
      BaseType;
  typedef StationaryContainerBasedDefault< ProblemImp, AnsatzSpaceImp, MatrixImp, VectorImp, TestSpaceImp > ThisType;
public:
  using typename BaseType::ProblemType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorType;

  StationaryContainerBasedDefault(const ProblemType& prblm,
                                  AnsatzSpaceType ansatz_sp,
                                  TestSpaceType test_sp,
                                  MatrixType system_mtrx,
                                  VectorType rhs_vec,
                                  VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {}

  StationaryContainerBasedDefault(const ProblemType& prblm,
                                  AnsatzSpaceType ansatz_sp,
                                  MatrixType system_mtrx,
                                  VectorType rhs_vec,
                                  VectorType dirichlet)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_(dirichlet)
    , has_dirichlet_shift_(true)
  {}

  StationaryContainerBasedDefault(const ProblemType& prblm,
                                  AnsatzSpaceType ansatz_sp,
                                  TestSpaceType test_sp,
                                  MatrixType system_mtrx,
                                  VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(test_sp)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {}

  StationaryContainerBasedDefault(const ProblemType& prblm,
                                  AnsatzSpaceType ansatz_sp,
                                  MatrixType system_mtrx,
                                  VectorType rhs_vec)
    : problem_(prblm)
    , ansatz_space_(ansatz_sp)
    , test_space_(ansatz_space_)
    , system_matrix_(system_mtrx)
    , rhs_vector_(rhs_vec)
    , dirichlet_shift_()
    , has_dirichlet_shift_(false)
  {}

  StationaryContainerBasedDefault(ThisType&& /*source*/) = default;

  /// \name Required by StationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return problem_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  /// \}
  /// \name Required by ContainerBasedStationaryDiscretizationInterface.
  /// \{

  const MatrixType& system_matrix() const
  {
    return system_matrix_;
  }

  const VectorType& rhs_vector() const
  {
    return rhs_vector_;
  }

  bool has_dirichlet_shift() const
  {
    return has_dirichlet_shift_;
  }

  const VectorType& dirichlet_shift() const
  {
    if (!has_dirichlet_shift_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call dirichlet_shift() if has_dirichlet_shift() is false!");
    return dirichlet_shift_;
  }

  /// \}

private:
  const ProblemType& problem_;
  const AnsatzSpaceType ansatz_space_;
  const TestSpaceType test_space_;
  const MatrixType system_matrix_;
  const VectorType rhs_vector_;
  const VectorType dirichlet_shift_;
  const bool has_dirichlet_shift_;
}; // class StationaryContainerBasedDefault


template< class ProblemImp, class FVSpaceImp >
class NonStationaryDefault
  : public NonStationaryDiscretizationInterface<
             internal::NonStationaryDefaultTraits< ProblemImp, FVSpaceImp > >
{
  typedef NonStationaryDiscretizationInterface
          < internal::NonStationaryDefaultTraits< ProblemImp, FVSpaceImp > >
      BaseType;
  typedef NonStationaryDefault< ProblemImp, FVSpaceImp > ThisType;
public:
  using typename BaseType::ProblemType;
  using typename BaseType::FVSpaceType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::VectorType;
  using typename BaseType::DiscreteFunctionType;

  NonStationaryDefault(const ProblemType& prblm,
                       const std::shared_ptr< const FVSpaceType > fv_space_ptr)
    : problem_(prblm)
    , fv_space_(fv_space_ptr)
  {}

  /// \name Required by NonStationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return problem_;
  }

  const FVSpaceType& fv_space() const
  {
    return *fv_space_;
  }

  using BaseType::solve;

  void solve(DiscreteSolutionType& solution, const bool is_linear) const
  {
#if HAVE_EIGEN
    try {
      // set dimensions
      static const size_t dimDomain = ProblemType::dimDomain;

      //get analytical flux, initial and boundary values
      typedef typename ProblemType::FluxType              AnalyticalFluxType;
      typedef typename ProblemType::RHSType               RHSType;
      typedef typename ProblemType::InitialValueType      InitialValueType;
      typedef typename ProblemType::BoundaryValueType     BoundaryValueType;
      typedef typename ProblemType::RangeFieldType        RangeFieldType;
      const std::shared_ptr< const AnalyticalFluxType > analytical_flux = problem_.flux();
      const std::shared_ptr< const InitialValueType > initial_values = problem_.initial_values();
      const std::shared_ptr< const BoundaryValueType > boundary_values = problem_.boundary_values();
      const std::shared_ptr< const RHSType > rhs = problem_.rhs();

      // allocate a discrete function for the concentration and another one to temporary store the update in each step
      typedef DiscreteFunction< FVSpaceType, Dune::Stuff::LA::CommonDenseVector< RangeFieldType > > FVFunctionType;
      FVFunctionType u(*fv_space_, "solution");

      //project initial values
      project(*initial_values, u);

      const double t_end = problem_.t_end();
      const double CFL = problem_.CFL();

      //calculate dx and choose t_end and initial dt
      Dune::Stuff::Grid::Dimensions< typename FVSpaceType::GridViewType > dimensions(fv_space_->grid_view());
      double dx = dimensions.entity_width.max();
      if (dimDomain == 2)
        dx /= std::sqrt(2);
      double dt = CFL*dx;

      // define operator types
      typedef typename Dune::GDT::Operators::AdvectionGodunov
          < AnalyticalFluxType, BoundaryValueType > OperatorType;
//      typedef typename Dune::Stuff::Functions::Constant< typename FVSpaceType::EntityType,
//                                                         DomainFieldType, dimDomain,
//                                                         RangeFieldType, 1, 1 >        ConstantFunctionType;
//      typedef typename Dune::GDT::Operators::AdvectionLaxFriedrichs
//          < AnalyticalFluxType, BoundaryValueType, ConstantFunctionType > OperatorType;
      typedef typename Dune::GDT::Operators::AdvectionRHS< RHSType > RHSOperatorType;

      // create right hand side operator
      RHSOperatorType rhs_operator(*rhs);

      //create advection operator
      OperatorType advection_operator(*analytical_flux, *boundary_values, is_linear);
//      const ConstantFunctionType dx_function(dx);
//      OperatorType advection_operator(*analytical_flux, *boundary_values, dx_function, dt, is_linear, false);

      typedef typename Dune::GDT::TimeStepper::ExplicitRungeKutta< OperatorType, FVFunctionType, double > OperatorTimeStepperType;
      OperatorTimeStepperType timestepper_op(advection_operator, u, -1.0);

      // do the time steps
      const size_t num_save_steps = 100;
      solution.clear();
      if (problem_.has_non_zero_rhs()) {
        // use fractional step method
        typedef typename Dune::GDT::TimeStepper::ExplicitRungeKutta< RHSOperatorType, FVFunctionType, double > RHSOperatorTimeStepperType;
        RHSOperatorTimeStepperType timestepper_rhs(rhs_operator, u);
        typedef typename Dune::GDT::TimeStepper::FractionalStep< OperatorTimeStepperType, RHSOperatorTimeStepperType > TimeStepperType;
        TimeStepperType timestepper(timestepper_op, timestepper_rhs);
        timestepper.solve(t_end, dt, num_save_steps, solution);
      } else {
        timestepper_op.solve(t_end, dt, num_save_steps, solution);
      }

    } catch (Dune::Exception& e) {
      std::cerr << "Dune reported: " << e.what() << std::endl;
      std::abort();
    }
#else //HAVE_EIGEN
    static_assert(AlwaysFalse< DiscreteSolutionType >::value, "You are missing eigen!");
#endif //HAVE_EIGEN
  }

  /// \}

private:
  const ProblemType& problem_;
  const std::shared_ptr< const FVSpaceType > fv_space_;
}; // class NonStationaryDefault



} // namespace Discretizations
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
