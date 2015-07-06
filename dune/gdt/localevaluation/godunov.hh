// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_EVALUATION_GODUNOV_HH
#define DUNE_GDT_EVALUATION_GODUNOV_HH

#include <tuple>
#include <memory>

# if HAVE_EIGEN
#   include <Eigen/Eigenvalues>
# endif

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace Godunov {


template< class LocalizableFunctionImp, size_t domainDim >
class Inner;

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
class Dirichlet;

template< class SourceFunctionImp >
class SourceEvaluation;

namespace internal {

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation.
 */
template< class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class InnerTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                                    LocalizableFunctionType;
  typedef Inner< LocalizableFunctionType, domainDim >               derived_type;
  typedef typename LocalizableFunctionType::EntityType              EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType         DomainFieldType;
  typedef typename LocalizableFunctionType::RangeFieldType          RangeFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType       LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >        LocalfunctionTupleType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = LocalizableFunctionType::dimRange;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity              FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType,
                                                RangeFieldType, dimRange,
                                                RangeFieldType, dimRange, dimDomain >   AnalyticalFluxType;

  typedef typename AnalyticalFluxType::RangeType                                        FluxRangeType;
  typedef typename AnalyticalFluxType::JacobianRangeType                                FluxJacobianRangeType;
  typedef typename Dune::Stuff::LocalfunctionSetInterface< EntityType,
                                                           DomainFieldType, dimDomain,
                                                           RangeFieldType, dimRange, 1 >::RangeType  RangeType;
  typedef typename Dune::Stuff::LA::EigenDenseMatrix< RangeFieldType >                  EigenMatrixType;
}; // class InnerTraits

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class DirichletTraits
   : public InnerTraits< LocalizableFunctionImp, domainDim >
{
  typedef InnerTraits< LocalizableFunctionImp, domainDim > BaseType;
public:
  typedef LocalizableFunctionImp                                            LocalizableFunctionType;
  typedef BoundaryValueFunctionImp                                          BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType             BoundaryValueLocalfunctionType;
  typedef Dirichlet< LocalizableFunctionType, BoundaryValueFunctionType, domainDim >   derived_type;
  using typename BaseType::LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType >,
                      std::shared_ptr< BoundaryValueLocalfunctionType > >   LocalfunctionTupleType;
}; // class DirichletTraits

template< class SourceFunctionImp >
class SourceEvaluationTraits
{
public:
  typedef SourceFunctionImp                                            SourceFunctionType;
  typedef SourceEvaluation< SourceFunctionType >                       derived_type;
  typedef typename SourceFunctionType::EntityType                      EntityType;
  typedef typename SourceFunctionType::LocalfunctionType               LocalfunctionType;
  typedef std::tuple< typename std::unique_ptr< LocalfunctionType > >  LocalfunctionTupleType;
  typedef typename SourceFunctionType::DomainFieldType                 DomainFieldType;
  static const size_t dimDomain = SourceFunctionType::dimDomain;
  static const size_t dimRange = SourceFunctionType::rangeDimRange;
};


} // namespace internal



#define PAPERFLUX 1

template< class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class Inner
  : public LocalEvaluation::Codim1Interface< internal::InnerTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::InnerTraits< LocalizableFunctionImp, domainDim >           Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit Inner(const AnalyticalFluxType& analytical_flux,
                 const LocalizableFunctionType& dx,
                 const double dt,
                 const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity));
  }

  size_t order(const LocalfunctionTupleType& /*localFunctionsEntity*/,
               const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*ansatzBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*ansatzBaseNeighbor*/) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template< class IntersectionType >
  void evaluate(const LocalfunctionTupleType& /*localFunctionsEntity*/,
                const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< RangeFieldType >& /*entityEntityRet*/,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< RangeFieldType >& entityNeighborRet,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborEntityRet*/) const
  {
    // get function values
    const RangeType u_i = ansatzBaseEntity.evaluate(intersection.geometryInInside().center())[0];
    const RangeType u_j = ansatzBaseNeighbor.evaluate(intersection.geometryInOutside().center())[0];
    // get flux value
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;
    // get unit outer normal
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
#ifndef NDEBUG
        ++num_zeros;
#endif NDEBUG
      }
      else
        DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
    }
    assert(num_zeros == dimDomain - 1);
    // calculate return vector
    RangeFieldType vol_intersection = intersection.geometry().volume();
    if (n_ij[coord] > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_[coord].mv(delta_u, negative_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        entityNeighborRet[0][kk] = (f_u_i[kk][coord] - negative_waves[kk]*n_ij[coord])*vol_intersection;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_[coord].mv(delta_u, positive_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        entityNeighborRet[0][kk] = (-f_u_i[kk][coord] - positive_waves[kk]*n_ij[coord])*vol_intersection;
    }
  } // void evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i,
                              const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector< EigenMatrixType >&& jacobian) const
  {
#if HAVE_EIGEN
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
          eigen_solver(jacobian[ii].backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      }
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
      jacobian_neg_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
      jacobian_pos_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Inner

template < class LocalizableFunctionImp, size_t dimDomain >
typename internal::InnerTraits< LocalizableFunctionImp, dimDomain >::FluxJacobianRangeType
Inner< LocalizableFunctionImp, dimDomain >::jacobian_neg_(0);

template < class LocalizableFunctionImp, size_t dimDomain >
typename internal::InnerTraits< LocalizableFunctionImp, dimDomain >::FluxJacobianRangeType
Inner< LocalizableFunctionImp, dimDomain >::jacobian_pos_(0);

template < class LocalizableFunctionImp, size_t dimDomain >
bool
Inner< LocalizableFunctionImp, dimDomain >::jacobians_constructed_(false);

template< class LocalizableFunctionImp >
class Inner< LocalizableFunctionImp, 1 >
  : public LocalEvaluation::Codim1Interface< internal::InnerTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::InnerTraits< LocalizableFunctionImp, 1 >        Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename DS::Functions::Affine< typename AnalyticalFluxType::EntityType,
                                          RangeFieldType, dimRange,
                                          RangeFieldType, dimRange, 1 >         AffineFunctionType;

  explicit Inner(const AnalyticalFluxType& analytical_flux,
                 const LocalizableFunctionType& dx,
                 const double dt,
                 const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity));
  }

  size_t order(const LocalfunctionTupleType& /*localFunctionsEntity*/,
               const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*ansatzBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*ansatzBaseNeighbor*/) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template< class IntersectionType >
  void evaluate(const LocalfunctionTupleType& /*localFunctionsEntity*/,
                const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< RangeFieldType >& /*entityEntityRet*/,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< RangeFieldType >& entityNeighborRet,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborEntityRet*/) const
  {
    // get function values
    const RangeType u_i = ansatzBaseEntity.evaluate(intersection.geometryInInside().center())[0];
    const RangeType u_j = ansatzBaseNeighbor.evaluate(intersection.geometryInOutside().center())[0];

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get unit outer normal
    const RangeFieldType n_ij = intersection.unitOuterNormal(localPoint)[0];
    // calculate return vector
#if PAPERFLUX
    //get flux values, FluxRangeType should be FieldVector< ..., dimRange >
    const FluxRangeType f_u_i_plus_f_u_j = is_linear_
                                           ? analytical_flux_.evaluate(u_i + u_j)
                                           : analytical_flux_.evaluate(u_i) + analytical_flux_.evaluate(u_j);
    RangeType waves(RangeFieldType(0));
    if (n_ij > 0) {
      // flux = 0.5*(f_u_i + f_u_j + |A|*(u_i-u_j))*n_ij
      jacobian_abs_function_->evaluate(u_i - u_j, waves);
      entityNeighborRet[0].axpy(RangeFieldType(0.5), f_u_i_plus_f_u_j + waves);
    } else {
      // flux = 0.5*(f_u_i + f_u_j - |A|*(u_i-u_j))*n_ij
      jacobian_abs_function_->evaluate(u_j - u_i, waves);
      entityNeighborRet[0].axpy(RangeFieldType(-0.5), f_u_i_plus_f_u_j + waves);
    }
#else
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
    if (n_ij > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_function_->evaluate(u_j - u_i, negative_waves);
      entityNeighborRet[0] = f_u_i + negative_waves;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_function_->evaluate(u_i - u_j, positive_waves);
      entityNeighborRet[0] = positive_waves - f_u_i;
    }
#endif
  } // void evaluate(...) const

private:
  void initialize_jacobians() const
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    EigenMatrixType jacobian_eigen(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian)));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i,
                              const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    EigenMatrixType jacobian_eigen = DSC::fromString< EigenMatrixType >(DSC::toString(jacobian));

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(EigenMatrixType&& jacobian) const
  {
#if HAVE_EIGEN
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    diag_jacobian_pos_tmp.scal(RangeFieldType(0));
    diag_jacobian_neg_tmp.scal(RangeFieldType(0));
    // create EigenSolver
    ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
        eigen_solver(jacobian.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
    assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
    for (size_t jj = 0; jj < dimRange; ++jj) {
      // assert this is real
      assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
      const RangeFieldType eigenvalue = eigenvalues[jj].real();
      if (eigenvalue < 0)
        diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
      else
        diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
      // set jacobian_neg_ and jacobian_pos_
      jacobian_neg_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
      jacobian_pos_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
      // jacobian_abs_ = jacobian_pos_ - jacobian_neg_;
      jacobian_abs_ = jacobian_neg_;
      jacobian_abs_ *= RangeFieldType(-1.0);
      jacobian_abs_ += jacobian_pos_;
# if PAPERFLUX
      jacobian_abs_function_ = DSC::make_unique< AffineFunctionType >(jacobian_abs_, RangeType(0), true);
# else
      jacobian_neg_function_ = DSC::make_unique< AffineFunctionType >(jacobian_neg_, RangeType(0), true);
      jacobian_pos_function_ = DSC::make_unique< AffineFunctionType >(jacobian_pos_, RangeType(0), true);
# endif
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static FluxJacobianRangeType jacobian_abs_;
  static std::unique_ptr< AffineFunctionType > jacobian_neg_function_;
  static std::unique_ptr< AffineFunctionType > jacobian_pos_function_;
  static std::unique_ptr< AffineFunctionType > jacobian_abs_function_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Inner< ..., 1 >

template < class LocalizableFunctionImp >
typename Inner< LocalizableFunctionImp, 1 >::FluxJacobianRangeType
Inner< LocalizableFunctionImp, 1 >::jacobian_neg_(0);

template < class LocalizableFunctionImp >
typename Inner< LocalizableFunctionImp, 1 >::FluxJacobianRangeType
Inner< LocalizableFunctionImp, 1 >::jacobian_pos_(0);

template < class LocalizableFunctionImp >
typename Inner< LocalizableFunctionImp, 1 >::FluxJacobianRangeType
Inner< LocalizableFunctionImp, 1 >::jacobian_abs_(0);

template < class LocalizableFunctionImp >
std::unique_ptr< typename Inner< LocalizableFunctionImp, 1 >::AffineFunctionType >
Inner< LocalizableFunctionImp, 1 >::jacobian_neg_function_;

template < class LocalizableFunctionImp >
std::unique_ptr< typename Inner< LocalizableFunctionImp, 1 >::AffineFunctionType >
Inner< LocalizableFunctionImp, 1 >::jacobian_pos_function_;

template < class LocalizableFunctionImp >
std::unique_ptr< typename Inner< LocalizableFunctionImp, 1 >::AffineFunctionType >
Inner< LocalizableFunctionImp, 1 >::jacobian_abs_function_;

template < class LocalizableFunctionImp >
bool
Inner< LocalizableFunctionImp, 1 >::jacobians_constructed_(false);

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class Dirichlet
  : public LocalEvaluation::Codim1Interface
                          < internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >, 2 >
{
public:
  typedef internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit Dirichlet(const AnalyticalFluxType& analytical_flux,
                     const LocalizableFunctionType& dx,
                     const double dt,
                     const BoundaryValueFunctionType& boundary_values,
                     const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity), boundary_values_.local_function(entity));
  }

  template< class R, unsigned long rT, unsigned long rCT, unsigned long rA, unsigned long rCA >
  size_t order(const LocalfunctionTupleType /*localFuncs*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBase*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBase*/) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
      const auto intersection_center_local = intersection.geometryInInside().center();
      const RangeType u_i = ansatzBase.evaluate(intersection_center_local)[0];
      const RangeType u_j = std::get< 1 >(localFunctions)->evaluate(intersection_center_local);
      // get flux values
      const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
      const FluxRangeType f_u_j = analytical_flux_.evaluate(u_j);

      if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
        reinitialize_jacobians(u_i, u_j);
      }
      // get jump at the intersection
      const RangeType delta_u = u_i - u_j;
      // get unit outer normal
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      // find direction of unit outer normal
      size_t coord = 0;
  #ifndef NDEBUG
      size_t num_zeros = 0;
  #endif NDEBUG
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
          coord = ii;
        else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
  #ifndef NDEBUG
          ++num_zeros;
  #endif NDEBUG
        }
        else
          DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
      }
      assert(num_zeros == dimDomain - 1);
      // calculate return vector
      RangeFieldType vol_intersection = intersection.geometry().volume();
      if (n_ij[coord] > 0) {
        RangeType negative_waves(RangeFieldType(0));
        jacobian_neg_[coord].mv(delta_u, negative_waves);
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[0][kk] = (f_u_i[kk][coord] - negative_waves[kk]*n_ij[coord])*vol_intersection;
      } else {
        RangeType positive_waves(RangeFieldType(0));
        jacobian_pos_[coord].mv(delta_u, positive_waves);
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[0][kk] = (-f_u_i[kk][coord] - positive_waves[kk]*n_ij[coord])*vol_intersection;
      }
    } // void evaluate(...) const

  private:
    void initialize_jacobians()
    {
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
      std::vector< EigenMatrixType > jacobian_eigen;
      for (size_t ii = 0; ii < dimDomain; ++ii)
        jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
      calculate_jacobians(std::move(jacobian_eigen));
      jacobians_constructed_ = true;
    } // void initialize_jacobians()

    void reinitialize_jacobians(const RangeType& u_i,
                                const RangeType& u_j) const
    {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      RangeType u_mean = u_i + u_j;
      u_mean *= RangeFieldType(0.5);
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
      std::vector< EigenMatrixType > jacobian_eigen;
      for (size_t ii = 0; ii < dimDomain; ++ii)
        jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
      // calculate jacobian_neg and jacobian_pos
      calculate_jacobians(std::move(jacobian_eigen));
    } // void reinitialize_jacobians(...)

    void calculate_jacobians(std::vector< EigenMatrixType >&& jacobian) const
    {
  #if HAVE_EIGEN
      EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
      EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        diag_jacobian_pos_tmp.scal(RangeFieldType(0));
        diag_jacobian_neg_tmp.scal(RangeFieldType(0));
        // create EigenSolver
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver(jacobian[ii].backend());
        assert(eigen_solver.info() == ::Eigen::Success);
        const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
        const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues[jj].real();
          if (eigenvalue < 0)
            diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
          else
            diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
        }
        const auto eigenvectors_inverse = eigenvectors.inverse();
        EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
        EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
        jacobian_neg_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
        jacobian_pos_[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
      }
  #else
      static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
  #endif
    } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const BoundaryValueFunctionType& boundary_values_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Dirichlet

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
typename internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >::FluxJacobianRangeType
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >::jacobian_neg_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
typename internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >::FluxJacobianRangeType
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >::jacobian_pos_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
bool
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >::jacobians_constructed_(false);

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >
  : public LocalEvaluation::Codim1Interface
                          < internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >, 2 >
{
public:
  typedef internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename DS::Functions::Affine< typename AnalyticalFluxType::EntityType,
                                          RangeFieldType, dimRange,
                                          RangeFieldType, dimRange, 1 >         AffineFunctionType;
  explicit Dirichlet(const AnalyticalFluxType& analytical_flux,
                     const LocalizableFunctionType& dx,
                     const double dt,
                     const BoundaryValueFunctionType& boundary_values,
                     const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity), boundary_values_.local_function(entity));
  }

  template< class R, unsigned long rT, unsigned long rCT, unsigned long rA, unsigned long rCA >
  size_t order(const LocalfunctionTupleType /*localFuncs*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBase*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBase*/) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
      const auto intersection_center_local = intersection.geometryInInside().center();
      const RangeType u_i = ansatzBase.evaluate(intersection_center_local)[0];
      const RangeType u_j = std::get< 1 >(localFunctions)->evaluate(intersection_center_local);

      if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
        reinitialize_jacobians(u_i, u_j);
      }
      // get unit outer normal
      const RangeFieldType n_ij = intersection.unitOuterNormal(localPoint)[0];
      // calculate return vector
  #if PAPERFLUX
      //get flux values, FluxRangeType should be FieldVector< ..., dimRange >
      const FluxRangeType f_u_i_plus_f_u_j = is_linear_
                                             ? analytical_flux_.evaluate(u_i + u_j)
                                             : analytical_flux_.evaluate(u_i) + analytical_flux_.evaluate(u_j);
      RangeType waves(RangeFieldType(0));
      if (n_ij > 0) {
        // flux = 0.5*(f_u_i + f_u_j + |A|*(u_i-u_j))*n_ij
        jacobian_abs_function_->evaluate(u_i - u_j, waves);
        ret[0].axpy(RangeFieldType(0.5), f_u_i_plus_f_u_j + waves);
      } else {
        // flux = 0.5*(f_u_i + f_u_j - |A|*(u_i-u_j))*n_ij
        jacobian_abs_function_->evaluate(u_j - u_i, waves);
        ret[0].axpy(RangeFieldType(-0.5), f_u_i_plus_f_u_j + waves);
      }
  #else
      const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
      if (n_ij > 0) {
        RangeType negative_waves(RangeFieldType(0));
        jacobian_neg_function_->evaluate(u_j - u_i, negative_waves);
        ret[0] = f_u_i + negative_waves;
      } else {
        RangeType positive_waves(RangeFieldType(0));
        jacobian_pos_function_->evaluate(u_i - u_j, positive_waves);
        ret[0] = positive_waves - f_u_i;
      }
  #endif
    }  // void evaluate(...) const

  private:
    void initialize_jacobians() const
    {
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
      EigenMatrixType jacobian_eigen(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian)));
      calculate_jacobians(jacobian_eigen);
      jacobians_constructed_ = true;
    } // void initialize_jacobians()

    void reinitialize_jacobians(const RangeType& u_i,
                                const RangeType& u_j) const
    {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      RangeType u_mean = u_i + u_j;
      u_mean *= RangeFieldType(0.5);
      const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
      EigenMatrixType jacobian_eigen = DSC::fromString< EigenMatrixType >(DSC::toString(jacobian));

      // calculate jacobian_neg and jacobian_pos
      calculate_jacobians(jacobian_eigen);
    } // void reinitialize_jacobians(...)

    void calculate_jacobians(EigenMatrixType& jacobian) const
    {
#if HAVE_EIGEN
      EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
      EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
          eigen_solver(jacobian.backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast< size_t >(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
        const auto eigenvectors_inverse = eigenvectors.inverse();
        EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
        EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
        // set jacobian_neg_ and jacobian_pos_
        jacobian_neg_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
        jacobian_pos_ = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
        // jacobian_abs_ = jacobian_pos_ - jacobian_neg_;
        jacobian_abs_ = jacobian_neg_;
        jacobian_abs_ *= RangeFieldType(-1.0);
        jacobian_abs_ += jacobian_pos_;
# if PAPERFLUX
      jacobian_abs_function_ = DSC::make_unique< AffineFunctionType >(jacobian_abs_, RangeType(0), true);
# else
      jacobian_neg_function_ = DSC::make_unique< AffineFunctionType >(jacobian_neg_, RangeType(0), true);
      jacobian_pos_function_ = DSC::make_unique< AffineFunctionType >(jacobian_pos_, RangeType(0), true);
# endif
      }
#else
      static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
    } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const BoundaryValueFunctionType& boundary_values_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static FluxJacobianRangeType jacobian_abs_;
  static std::unique_ptr< AffineFunctionType > jacobian_neg_function_;
  static std::unique_ptr< AffineFunctionType > jacobian_pos_function_;
  static std::unique_ptr< AffineFunctionType > jacobian_abs_function_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Dirichlet< ..., 1 >

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_neg_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_pos_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::FluxJacobianRangeType
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_abs_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
std::unique_ptr< typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType >
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_neg_function_;

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
std::unique_ptr< typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType >
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_pos_function_;

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
std::unique_ptr< typename Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::AffineFunctionType >
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobian_abs_function_;

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
bool
Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >::jacobians_constructed_(false);


template< class SourceFunctionImp >
class SourceEvaluation
  : public Codim0Interface< internal::SourceEvaluationTraits< SourceFunctionImp >, 1 >
{
public:
  typedef typename internal::SourceEvaluationTraits< SourceFunctionImp > Traits;
  typedef typename Traits::SourceFunctionType                            SourceFunctionType;
  typedef typename Traits::EntityType                                    EntityType;
  typedef typename Traits::LocalfunctionTupleType                        LocalfunctionTupleType;
  typedef typename Traits::DomainFieldType                               DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit SourceEvaluation(const SourceFunctionType& source_function)
    : source_function_(source_function)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(source_function_.local_global_function(entity));
  }

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class R, size_t r, size_t rC >
  size_t order(const LocalfunctionTupleType& /*localFunctions_in*/,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& /*testBase*/)
  const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeCols of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class R, size_t rC >
  void evaluate(const LocalfunctionTupleType& local_source_function,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, rC >& entityAverage,
                const Dune::FieldVector< DomainFieldType, dimDomain >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
#if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) //EXADUNE
    ret = std::get< 0 >(local_source_function)->evaluate(localPoint, entityAverage.evaluate(localPoint)[0]);
#else
    const auto fieldvector_ret = std::get< 0 >(local_source_function)->evaluate(localPoint,
                                                                                entityAverage.evaluate(localPoint)[0]);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = fieldvector_ret[ii];
#endif
  }

private:
  const SourceFunctionType& source_function_;
}; // class Codim0Interface< Traits, 1 >


} // namespace Godunov
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_GODUNOV_HH
