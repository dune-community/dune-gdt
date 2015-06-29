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

#include <dune/gdt/operators/advection.hh>
#include <dune/gdt/timestepper/rungekutta.hh>

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
  typedef typename Dune::Stuff::LA::CommonDenseVector< RangeFieldType >   StationaryVectorType;
  typedef DiscreteFunction< FVSpaceType, StationaryVectorType >           DiscreteFunctionType;
  typedef std::vector< std::pair< double, StationaryVectorType > >        VectorType;
}; // class StationaryContainerBasedDefaultTraits


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
  using typename BaseType::VectorType;
  using typename BaseType::StationaryVectorType;
  using typename BaseType::DiscreteFunctionType;

  NonStationaryDefault(const ProblemType& prblm,
                       const FVSpaceType fvspace)
    : problem_(prblm)
    , fv_space_(fvspace)
  {}

  NonStationaryDefault(ThisType&& /*source*/) = default;

  /// \name Required by NonStationaryDiscretizationInterface.
  /// \{

  const ProblemType& problem() const
  {
    return problem_;
  }

  const FVSpaceType& fv_space() const
  {
    return fv_space_;
  }

  using BaseType::solve;

  void solve(VectorType& solution, const bool is_linear/* = false*/) const
  {
    try {
      static const size_t dimDomain = ProblemType::dimDomain;
      static const size_t dimRange = ProblemType::dimRange;

      //get analytical flux and initial values
      typedef typename ProblemType::FluxType            AnalyticalFluxType;
      typedef typename ProblemType::SourceType          SourceType;
      typedef typename ProblemType::FunctionType        FunctionType;
      typedef typename ProblemType::BoundaryValueType   BoundaryValueType;
      typedef typename FunctionType::DomainFieldType    DomainFieldType;
      typedef typename ProblemType::RangeFieldType      RangeFieldType;
      const std::shared_ptr< const FunctionType > initial_values = problem_.initial_values();
      const std::shared_ptr< const AnalyticalFluxType > analytical_flux = problem_.flux();
      const std::shared_ptr< const BoundaryValueType > boundary_values = problem_.boundary_values();
      const std::shared_ptr< const SourceType > source = problem_.source();

      // allocate a discrete function for the concentration and another one to temporary store the update in each step
      DiscreteFunctionType u(fv_space_, "solution");

      //project initial values
      project(*initial_values, u);

      const double t_end = problem_.t_end();
      const double ratio_dt_dx = problem_.ratio_dt_dx();

      //calculate h and then dt from the fixed CFL ratio
      Dune::Stuff::Grid::Dimensions< typename FVSpaceType::GridViewType > dimensions(fv_space_.grid_view());
      const double dx = dimensions.entity_width.max();
      const double dt = ratio_dt_dx*dx;
      typedef typename Dune::Stuff::Functions::Constant< typename FVSpaceType::EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 > ConstantFunctionType;
      ConstantFunctionType dx_function(dx);

      //create operator
      typedef typename Dune::GDT::Operators::AdvectionGodunov< AnalyticalFluxType, ConstantFunctionType, BoundaryValueType, FVSpaceType/*, Operators::SlopeLimiters::mc*/ > OperatorType;
      OperatorType advection_operator(*analytical_flux, dx_function, dt, *boundary_values, fv_space_, is_linear);

      //create butcher_array
      // forward euler
      Dune::DynamicMatrix< RangeFieldType > A(DSC::fromString< Dune::DynamicMatrix< RangeFieldType >  >("[0]"));
      Dune::DynamicVector< RangeFieldType > b(DSC::fromString< Dune::DynamicVector< RangeFieldType >  >("[1]"));
      // generic second order, x = 1 (see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)
      //    Dune::DynamicMatrix< RangeFieldType > A(DSC::fromString< Dune::DynamicMatrix< RangeFieldType >  >("[0 0; 1 0]"));
      //    Dune::DynamicVector< RangeFieldType > b(DSC::fromString< Dune::DynamicVector< RangeFieldType >  >("[0.5 0.5]"));
      // classic fourth order RK
      //    Dune::DynamicMatrix< RangeFieldType > A(DSC::fromString< Dune::DynamicMatrix< RangeFieldType >  >("[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]"));
      //    Dune::DynamicVector< RangeFieldType > b(DSC::fromString< Dune::DynamicVector< RangeFieldType >  >("[" + DSC::toString(1.0/6.0) + " " + DSC::toString(1.0/3.0) + " " + DSC::toString(1.0/3.0) + " " + DSC::toString(1.0/6.0) + "]"));

      //create timestepper
      Dune::GDT::TimeStepper::RungeKutta< OperatorType, DiscreteFunctionType, SourceType > timestepper(advection_operator, u, *source, dx, A, b);

      // now do the time steps
      std::vector< std::pair< double, DiscreteFunctionType > > solution_as_discrete_function;

      const double saveInterval = t_end/500 > dt ? t_end/500 : dt;
      timestepper.solve(t_end, dt, saveInterval, false, false, true, solution_as_discrete_function);
      solution.clear();
      const size_t num_time_steps = solution_as_discrete_function.size();
      for (size_t ii = 0; ii < num_time_steps; ++ii) {
        StationaryVectorType stationary_vector(solution_as_discrete_function[ii].second.vector());
        solution.emplace_back(std::make_pair(solution_as_discrete_function[ii].first, stationary_vector));
      }
    } catch (Dune::Exception& e) {
      std::cerr << "Dune reported: " << e.what() << std::endl;
      std::abort();
    }
  }

  /// \}

private:
  const ProblemType& problem_;
  const FVSpaceType fv_space_;
}; // class NonStationaryDefault



} // namespace Discretizations
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETIZATIONS_DEFAULT_HH
