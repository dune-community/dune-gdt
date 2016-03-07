// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_FV_HH
#define DUNE_GDT_OPERATORS_FV_HH

#include <memory>
#include <type_traits>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/information.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localfluxes/interfaces.hh>
#include <dune/gdt/localfluxes/godunov.hh>
#include <dune/gdt/localfluxes/laxfriedrichs.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localoperator/fv.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/default.hh>

#include <dune/gdt/playground/spaces/dg/pdelabproduct.hh>

#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forwards
template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
class AdvectionLaxFriedrichs;

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter >
class AdvectionGodunov;

template< class RHSEvaluationImp >
class AdvectionRHS;


namespace internal {


//TODO: add static assert once type of BoundaryValueFunctionImp is decided
template< class AnalyticalFluxImp, class BoundaryValueFunctionImp >
class AdvectionTraitsBase
{
  static_assert(is_analytical_flux< AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
//  static_assert(Stuff::is_???< BoundaryValueFunctionImp >::value,
//                "BoundaryValueFunctionImp has to be derived from ???!");
public:
  typedef AnalyticalFluxImp                                                                     AnalyticalFluxType;
  typedef BoundaryValueFunctionImp                                                              BoundaryValueFunctionType;
  typedef typename AnalyticalFluxType::DomainFieldType                                          FieldType;
}; // class AdvectionTraitsBase


template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
class AdvectionLaxFriedrichsTraits
    : public AdvectionTraitsBase< AnalyticalFluxImp, BoundaryValueFunctionImp >
{
  static_assert(Stuff::is_localizable_function< LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef AdvectionLaxFriedrichs< AnalyticalFluxImp,
                                  BoundaryValueFunctionImp,
                                  LocalizableFunctionImp>                derived_type;

}; // class AdvectionLaxFriedrichsTraits

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter >
class AdvectionGodunovTraits
    : public AdvectionTraitsBase< AnalyticalFluxImp, BoundaryValueFunctionImp >
{
public:
  typedef AdvectionGodunov< AnalyticalFluxImp, BoundaryValueFunctionImp, slope_limiter > derived_type;
}; // class AdvectionGodunovTraits


template< class RHSEvaluationImp >
class AdvectionRHSTraits
{
public:
  typedef AdvectionRHS< RHSEvaluationImp >                                                        derived_type;
  typedef RHSEvaluationImp                                                                        RHSEvaluationType;
  typedef typename RHSEvaluationImp::DomainFieldType                                              FieldType;
}; // class AdvectionRHSTraits


} // namespace internal


template< class AnalyticalFluxImp,
          class NumericalCouplingFluxImp,
          class NumericalBoundaryFluxImp,
          class BoundaryValueFunctionImp,
          class SourceImp,
          class RangeImp >
class AdvectionLocalizableDefault
  : public Dune::GDT::LocalizableOperatorDefault< typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp >
{
  typedef Dune::GDT::LocalizableOperatorDefault< typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp > BaseType;

  static_assert(is_analytical_flux< AnalyticalFluxImp >::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface!");
  static_assert(is_numerical_coupling_flux< NumericalCouplingFluxImp >::value,
                "NumericalCouplingFluxImp has to be derived from NumericalCouplingFluxInterface!");
  static_assert(is_numerical_boundary_flux< NumericalBoundaryFluxImp >::value,
                "NumericalBoundaryFluxImp has to be derived from NumericalBoundaryFluxInterface!");
//  static_assert(std::is_base_of< ???, BoundaryValueFunctionImp >::value,
//                "BoundaryValueFunctionImp has to be derived from ???!");
  static_assert(is_discrete_function< SourceImp >::value, "SourceImp has to be derived from DiscreteFunction!");
  static_assert(is_discrete_function< RangeImp >::value, "RangeImp has to be derived from DiscreteFunction!");
public:
  typedef AnalyticalFluxImp                                            AnalyticalFluxType;
  typedef NumericalCouplingFluxImp                                     NumericalCouplingFluxType;
  typedef NumericalBoundaryFluxImp                                     NumericalBoundaryFluxType;
  typedef BoundaryValueFunctionImp                                     BoundaryValueFunctionType;
  typedef SourceImp                                                    SourceType;
  typedef RangeImp                                                     RangeType;
  typedef typename SourceType::RangeFieldType                          RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType                  GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalCouplingFVOperator< NumericalCouplingFluxType >   LocalCouplingOperatorType;
  typedef typename Dune::GDT::LocalBoundaryFVOperator< NumericalBoundaryFluxType >   LocalBoundaryOperatorType;

  template< class... LocalOperatorArgTypes >
  AdvectionLocalizableDefault(const AnalyticalFluxType& analytical_flux,
                              const BoundaryValueFunctionType& boundary_values,
                              const SourceType& source,
                              RangeType& range,
                              LocalOperatorArgTypes&&... local_operator_args)
    : BaseType(range.space().grid_view(), source, range)
    , local_operator_(analytical_flux, std::forward< LocalOperatorArgTypes >(local_operator_args)...)
    , local_boundary_operator_(analytical_flux, boundary_values, std::forward< LocalOperatorArgTypes >(local_operator_args)...)
  {
    this->add(local_operator_, new DSG::ApplyOn::InnerIntersectionsPrimally< GridViewType >());
    this->add(local_operator_, new DSG::ApplyOn::PeriodicIntersectionsPrimally< GridViewType >());
    this->add(local_boundary_operator_, new DSG::ApplyOn::NonPeriodicBoundaryIntersections< GridViewType >());
  }

private:
  const LocalCouplingOperatorType local_operator_;
  const LocalBoundaryOperatorType local_boundary_operator_;
}; // class AdvectionLocalizableDefault


template< class SourceImp, class RangeImp, class BoundaryValueFunctionImp, class MatrixImp, SlopeLimiters slope_limiter >
class LinearReconstructionLocalizable
  : public Dune::GDT::LocalizableOperatorDefault< typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp >
{
  typedef Dune::GDT::LocalizableOperatorDefault< typename RangeImp::SpaceType::GridViewType, SourceImp, RangeImp > BaseType;
  typedef LinearReconstructionLocalizable<SourceImp, RangeImp, BoundaryValueFunctionImp, MatrixImp, slope_limiter > ThisType;
public:
  typedef SourceImp                                                    SourceType;
  typedef RangeImp                                                     RangeType;
  typedef BoundaryValueFunctionImp                                     BoundaryValueFunctionType;
  typedef MatrixImp                                                    MatrixType;
  typedef typename SourceType::RangeFieldType                          RangeFieldType;
  typedef typename RangeType::SpaceType::GridViewType                  GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef typename Dune::GDT::LocalReconstructionFVOperator< MatrixType, BoundaryValueFunctionType, slope_limiter > LocalOperatorType;

  LinearReconstructionLocalizable(const SourceType& source,
                                  RangeType& range,
                                  const MatrixType& eigenvectors,
                                  const MatrixType& eigenvectors_inverse,
                                  const BoundaryValueFunctionType& boundary_values)
    : BaseType(range.space().grid_view(), source, range)
    , local_operator_(eigenvectors, eigenvectors_inverse, boundary_values)
    , source_(source)
    , range_(range)
  {
    this->add(local_operator_);
  }

private:
  const LocalOperatorType local_operator_;
  const SourceType& source_;
  RangeType& range_;
}; // class LinearReconstructionLocalizable


template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
class AdvectionLaxFriedrichs
  : public Dune::GDT::OperatorInterface< internal::AdvectionLaxFriedrichsTraits<  AnalyticalFluxImp,
                                                                                  BoundaryValueFunctionImp,
                                                                                  LocalizableFunctionImp > >
{
public:
  typedef internal::AdvectionLaxFriedrichsTraits< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp > Traits;
  typedef typename Traits::AnalyticalFluxType      AnalyticalFluxType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  typedef typename Dune::GDT::LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxType, LocalizableFunctionType, dimDomain >   NumericalCouplingFluxType;
  typedef typename Dune::GDT::LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxType, typename BoundaryValueFunctionType::TimeIndependentFunctionType, LocalizableFunctionType, dimDomain > NumericalBoundaryFluxType;

  AdvectionLaxFriedrichs(const AnalyticalFluxType& analytical_flux,
                         const BoundaryValueFunctionType& boundary_values,
                         const LocalizableFunctionType& dx,
                         const double dt,
                         const bool is_linear = false,
                         const bool use_local = false,
                         const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_local_(use_local)
    , entity_geometries_equal_(entity_geometries_equal)
  {}

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    const auto current_boundary_values = boundary_values_.evaluate_at_time(time);
    AdvectionLocalizableDefault< AnalyticalFluxType,
                                 NumericalCouplingFluxType,
                                 NumericalBoundaryFluxType,
                                 typename BoundaryValueFunctionType::ExpressionFunctionType,
                                 SourceType,
                                 RangeType
        > localizable_operator(analytical_flux_, *current_boundary_values, source, range, dx_, dt_, is_linear_, use_local_, entity_geometries_equal_);
    localizable_operator.apply();
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_local_;
  const bool entity_geometries_equal_;
}; // class AdvectionLaxFriedrichs

// TODO: remove eigen dependency of GodunovNumericalCouplingFlux/GodunovNumericalBoundaryFlux
#if HAVE_EIGEN

namespace internal {
  template< bool dimDomain_equals_1 >
  struct InitializerChooser
  {};
}

// TODO: 0 boundary by default, so no need to specify boundary conditions for periodic grid views
template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, SlopeLimiters slope_limiter = SlopeLimiters::minmod >
class AdvectionGodunov
  : public Dune::GDT::OperatorInterface< internal::AdvectionGodunovTraits<  AnalyticalFluxImp,
                                                                            BoundaryValueFunctionImp,
                                                                            slope_limiter > >
{
public:
  typedef internal::AdvectionGodunovTraits< AnalyticalFluxImp,
                                            BoundaryValueFunctionImp,
                                            slope_limiter >          Traits;
  typedef typename Traits::AnalyticalFluxType                     AnalyticalFluxType;
  typedef typename Traits::BoundaryValueFunctionType              BoundaryValueFunctionType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
  static const size_t dimRangeCols = AnalyticalFluxType::dimRangeCols;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
private:
  typedef typename Dune::Stuff::LA::EigenDenseMatrix< RangeFieldType > EigenMatrixType;
  typedef typename Dune::Stuff::Common::FieldMatrix< RangeFieldType, dimRange, dimRange > MatrixType;
public:
  typedef typename Dune::GDT::GodunovNumericalCouplingFlux< AnalyticalFluxType, dimDomain >   NumericalCouplingFluxType;
  typedef typename Dune::GDT::GodunovNumericalBoundaryFlux< AnalyticalFluxType, typename BoundaryValueFunctionType::TimeIndependentFunctionType, dimDomain > NumericalBoundaryFluxType;

  AdvectionGodunov(const AnalyticalFluxType& analytical_flux,
                   const BoundaryValueFunctionType& boundary_values,
                   const bool is_linear = false,
                   const bool use_linear_reconstruction = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
    , use_linear_reconstruction_(use_linear_reconstruction)
  {
     initialize(internal::InitializerChooser< dimDomain == 1 >());
  }

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double time = 0.0) const
  {
    const auto current_boundary_values = boundary_values_.evaluate_at_time(time);
    if (use_linear_reconstruction_)
    {
      typedef Spaces::DG::PdelabBasedProduct< typename SourceType::SpaceType::GridViewType,
                                              1, //polOrder
                                              RangeFieldType,
                                              dimRange,
                                              dimRangeCols >               DGSpaceType;
      typedef DiscreteFunction< DGSpaceType, typename SourceType::VectorType >  ReconstructedDiscreteFunctionType;
      const auto dg_space_ = DSC::make_unique< DGSpaceType >(range.space().grid_view());
      const auto reconstruction = DSC::make_unique< ReconstructedDiscreteFunctionType >(*dg_space_, "reconstructed");
      LinearReconstructionLocalizable
          < SourceType, ReconstructedDiscreteFunctionType, typename BoundaryValueFunctionType::ExpressionFunctionType, MatrixType, slope_limiter>
                    reconstruction_operator(source,
                                            *reconstruction,
                                            eigenvectors_,
                                            eigenvectors_inverse_,
                                            *current_boundary_values);
      reconstruction_operator.apply();
      AdvectionLocalizableDefault< AnalyticalFluxType,
          NumericalCouplingFluxType,
          NumericalBoundaryFluxType,
          typename BoundaryValueFunctionType::ExpressionFunctionType,
          ReconstructedDiscreteFunctionType,
          RangeType
          > localizable_operator(analytical_flux_, *current_boundary_values, *reconstruction, range, is_linear_);
      localizable_operator.apply();
    } else {
      AdvectionLocalizableDefault< AnalyticalFluxType,
          NumericalCouplingFluxType,
          NumericalBoundaryFluxType,
          typename BoundaryValueFunctionType::ExpressionFunctionType,
          SourceType,
          RangeType
          > localizable_operator(analytical_flux_, *current_boundary_values, source, range, is_linear_);
      localizable_operator.apply();
    }
  }

private:
  void initialize(const internal::InitializerChooser< false >&)
  {
    if (use_linear_reconstruction_) {
      assert(false && "Linear reconstruction is only implemented in 1D!");
    }
  }

  void initialize(const internal::InitializerChooser< true >&)
  {
    if (use_linear_reconstruction_) {
      assert(is_linear_ && "Linear reconstruction is only implemented for linear analytical fluxes!");
      // calculate matrix of eigenvectors of A, where A is the jacobian of the linear analytical flux, i.e. u_t + A*u_x = 0.
      // As the analytical flux is linear, the jacobian A is constant, so it is enough to evaluate at 0.
      ::Eigen::EigenSolver< typename EigenMatrixType::BackendType > eigen_solver(DSC::from_string< EigenMatrixType >(DSC::to_string(analytical_flux_.jacobian(typename AnalyticalFluxType::RangeType(0)))).backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigen_eigenvectors = eigen_solver.eigenvectors();
#ifndef NDEBUG
      for (size_t ii = 0; ii < dimRange; ++ii)
        for (size_t jj = 0; jj < dimRange; ++jj)
          assert(eigen_eigenvectors(ii,jj).imag() < 1e-15);
#endif
      eigenvectors_ = DSC::from_string< MatrixType >(DSC::to_string(EigenMatrixType(eigen_eigenvectors.real())));
      eigenvectors_inverse_ = DSC::from_string< MatrixType >(DSC::to_string(EigenMatrixType(eigen_eigenvectors.inverse().real())));
    }
  }


  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType&  boundary_values_;
  const bool is_linear_;
  const bool use_linear_reconstruction_;
  MatrixType eigenvectors_;
  MatrixType eigenvectors_inverse_;
}; // class AdvectionGodunov

#else // HAVE_EIGEN

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp >
class AdvectionGodunov
{
  static_assert(AlwaysFalse< AnalyticalFluxImp >::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN


template< class RHSEvaluationImp >
class AdvectionRHS
  : public Dune::GDT::OperatorInterface< internal::AdvectionRHSTraits< RHSEvaluationImp > >
{
  static_assert(is_rhs_evaluation< RHSEvaluationImp >::value, "RHSEvaluationImp has to be derived from RHSInterface!");
public:
  typedef internal::AdvectionRHSTraits< RHSEvaluationImp >          Traits;
  typedef typename Traits::RHSEvaluationType                        RHSEvaluationType;
  typedef LocalRHSFVOperator< RHSEvaluationType >                   LocalOperatorType;

  AdvectionRHS(const RHSEvaluationType& rhs_evaluation)
    : local_operator_(rhs_evaluation)
  {}

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range, const double /*time*/ = 0.0) const
  {
    LocalizableOperatorDefault< typename RangeType::SpaceType::GridViewType, SourceType, RangeType > localizable_operator(range.space().grid_view(), source, range);
    localizable_operator.add(local_operator_);
    localizable_operator.apply();
  }

private:
  const LocalOperatorType local_operator_;
}; // class AdvectionRHS


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_HH
