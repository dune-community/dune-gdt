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
#include <dune/common/fmatrixev.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace Godunov {


template< class LocalizableFunctionImp >
class Inner;

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet;

namespace internal {

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation.
 */
template< class LocalizableFunctionImp >
class InnerTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                                    LocalizableFunctionType;
  typedef Inner< LocalizableFunctionType >                          derived_type;
  typedef typename LocalizableFunctionType::EntityType              EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType         DomainFieldType;
  typedef typename LocalizableFunctionType::RangeFieldType          RangeFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType       LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >        LocalfunctionTupleType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
  static const size_t dimRange = LocalizableFunctionType::dimRange;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity              FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType,
                                                RangeFieldType, dimRange,
                                                RangeFieldType, dimRange, dimDomain >   AnalyticalFluxType;

  typedef typename AnalyticalFluxType::RangeType                                        FluxRangeType;
  typedef typename Dune::FieldVector< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange >, dimDomain > FluxJacobianRangeType;
  typedef typename Dune::Stuff::LocalfunctionSetInterface< EntityType,
                                                           DomainFieldType, dimDomain,
                                                           RangeFieldType, dimRange, 1 >::RangeType  RangeType;
  typedef typename Dune::Stuff::LA::EigenDenseMatrix< RangeFieldType >                  EigenMatrixType;
}; // class InnerTraits

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class DirichletTraits
   : public InnerTraits< LocalizableFunctionImp >
{
  typedef InnerTraits< LocalizableFunctionImp > BaseType;
public:
  typedef LocalizableFunctionImp                                            LocalizableFunctionType;
  typedef BoundaryValueFunctionImp                                          BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType             BoundaryValueLocalfunctionType;
  typedef Dirichlet< LocalizableFunctionType, BoundaryValueFunctionType >   derived_type;
  using typename BaseType::LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType >,
                      std::shared_ptr< BoundaryValueLocalfunctionType > >   LocalfunctionTupleType;
}; // class DirichletTraits


} // namespace internal


template< class LocalizableFunctionImp >
class Inner
  : public LocalEvaluation::Codim1Interface< internal::InnerTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::InnerTraits< LocalizableFunctionImp >           Traits;
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

  explicit Inner(const AnalyticalFluxType& analytical_flux,
                 const LocalizableFunctionType& ratio_dt_dx,
                 const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
    jacobians_constructed_ = true;
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(ratio_dt_dx_.local_function(entity));
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
    const EntityType& entity = ansatzBaseEntity.entity();
    const EntityType& neighbor = ansatzBaseNeighbor.entity();
    const std::vector< RangeType > u_i
                                 = ansatzBaseEntity.evaluate(entity.geometry().local(entity.geometry().center()));
    const std::vector< RangeType > u_j
                                 = ansatzBaseNeighbor.evaluate(neighbor.geometry().local(neighbor.geometry().center()));
    FluxJacobianRangeType jacobian_pos = jacobian_pos_;
    FluxJacobianRangeType jacobian_neg = jacobian_neg_;
    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i[0], u_j[0], jacobian_neg, jacobian_pos);
    }
    const RangeType delta_u = u_i[0] - u_j[0];
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    // find direction of unit outer normal
    size_t coord = 0;
    size_t num_zeros = 0;
    for (size_t ii = 0; ii < n_ij.size(); ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
        ++num_zeros;
      else
        DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
    }
    assert(num_zeros == dimRange - 1);
    // calculate return vector
    RangeType negative_waves(RangeFieldType(0));
    RangeType positive_waves(RangeFieldType(0));
    jacobian_neg[coord].mv(delta_u, negative_waves);
    jacobian_pos[coord].mv(delta_u, positive_waves);
//    if (n_ij > 0) {
//      for (size_t kk = 0; kk < dimRange; ++kk)
//        entityNeighborRet[kk][0] = (f_u_i[kk] + f_u_j[kk] + (positive_waves[kk] - negative_waves[kk]))*n_ij*0.5;
//    } else {
//      for (size_t kk = 0; kk < dimRange; ++kk)
//        entityNeighborRet[kk][0] = (f_u_i[kk] +  f_u_j[kk] - (positive_waves[kk] - negative_waves[kk]))*n_ij*0.5;
//    }
    if (n_ij[coord] > 0) {
      for (size_t kk = 0; kk < dimRange; ++kk)
        entityNeighborRet[kk][0] = -negative_waves[kk]*n_ij[coord];
    } else {
      for (size_t kk = 0; kk < dimRange; ++kk)
        entityNeighborRet[kk][0] = -positive_waves[kk]*n_ij[coord];
    }
  } // void evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    calculate_jacobians(jacobian_eigen, jacobian_neg_, jacobian_pos_);
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType u_i,
                              const RangeType u_j,
                              FluxJacobianRangeType& jacobian_neg,
                              FluxJacobianRangeType& jacobian_pos) const
  {
    // clear jacobian_pos and jacobian_neg
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      jacobian_neg[ii] *= RangeFieldType(0);
      jacobian_pos[ii] *= RangeFieldType(0);
    }
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i;
    for (size_t ii = 0; ii < u_mean.size(); ++ii) {
      u_mean[ii] += u_j[ii];
      u_mean[ii] *= 0.5;
    }
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    //    std::cout << DSC::toString(jacobian) << std::endl;
    //    std::cout << DSC::toString(u_mean) << std::endl;
    //    std::cout << DSC::toString(u_i) << "   " << DSC::toString(u_j) << std::endl;

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(jacobian_eigen, jacobian_neg, jacobian_pos);
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector< EigenMatrixType > jacobian,
                           FluxJacobianRangeType& jacobian_neg,
                           FluxJacobianRangeType& jacobian_pos) const
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
      jacobian_neg[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
      jacobian_pos[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Inner

template < class LocalizableFunctionImp >
typename internal::InnerTraits< LocalizableFunctionImp >::FluxJacobianRangeType Inner< LocalizableFunctionImp >::jacobian_neg_(0);

template < class LocalizableFunctionImp >
typename internal::InnerTraits< LocalizableFunctionImp >::FluxJacobianRangeType Inner< LocalizableFunctionImp >::jacobian_pos_(0);

template < class LocalizableFunctionImp >
bool Inner< LocalizableFunctionImp >::jacobians_constructed_(false);

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet
  : public LocalEvaluation::Codim1Interface
                          < internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >, 2 >
{
public:
  typedef internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >  Traits;
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
                     const LocalizableFunctionType& ratio_dt_dx,
                     const BoundaryValueFunctionType& boundary_values,
                     const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (!jacobians_constructed_)
      initialize_jacobians();
    jacobians_constructed_ = true;
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(ratio_dt_dx_.local_function(entity), boundary_values_.local_function(entity));
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
      const EntityType& entity = ansatzBase.entity();
      const std::vector< RangeType > u_i = ansatzBase.evaluate(entity.geometry().local(entity.geometry().center()));
      const auto local_center_intersection = entity.geometry().local(intersection.geometry().center());
      const auto& u_j = std::get< 1 >(localFunctions)->evaluate(local_center_intersection);
      FluxJacobianRangeType jacobian_pos = jacobian_pos_;
      FluxJacobianRangeType jacobian_neg = jacobian_neg_;
      if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
        reinitialize_jacobians(u_i[0], u_j, jacobian_neg, jacobian_pos);
      }
      const RangeType delta_u = u_i[0] - u_j;
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      // find direction of unit outer normal
      size_t coord = 0;
      size_t num_zeros = 0;
      for (size_t ii = 0; ii < n_ij.size(); ++ii) {
        if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
          coord = ii;
        else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
          ++num_zeros;
        else
          DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
      }
      assert(num_zeros == dimRange - 1);
      // calculate return vector
      RangeType negative_waves(RangeFieldType(0));
      RangeType positive_waves(RangeFieldType(0));
      jacobian_neg[coord].mv(delta_u, negative_waves);
      jacobian_pos[coord].mv(delta_u, positive_waves);
  //    if (n_ij > 0) {
  //      for (size_t kk = 0; kk < dimRange; ++kk)
  //        entityNeighborRet[kk][0] = (f_u_i[kk] + f_u_j[kk] + (positive_waves[kk] - negative_waves[kk]))*n_ij*0.5;
  //    } else {
  //      for (size_t kk = 0; kk < dimRange; ++kk)
  //        entityNeighborRet[kk][0] = (f_u_i[kk] +  f_u_j[kk] - (positive_waves[kk] - negative_waves[kk]))*n_ij*0.5;
  //    }
      if (n_ij[coord] > 0) {
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[kk][0] = -negative_waves[kk]*n_ij[coord];
      } else {
        for (size_t kk = 0; kk < dimRange; ++kk)
          ret[kk][0] = -positive_waves[kk]*n_ij[coord];
      }
  } // void evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    calculate_jacobians(jacobian_eigen, jacobian_neg_, jacobian_pos_);
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType u_i,
                              const RangeType u_j,
                              FluxJacobianRangeType& jacobian_neg,
                              FluxJacobianRangeType& jacobian_pos) const
  {
    // clear jacobian_pos and jacobian_neg
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      jacobian_neg[ii] *= RangeFieldType(0);
      jacobian_pos[ii] *= RangeFieldType(0);
    }
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i;
    for (size_t ii = 0; ii < u_mean.size(); ++ii) {
      u_mean[ii] += u_j[ii];
      u_mean[ii] *= 0.5;
    }
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector< EigenMatrixType > jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::fromString< EigenMatrixType >(DSC::toString(jacobian[ii])));
    //    std::cout << DSC::toString(jacobian) << std::endl;
    //    std::cout << DSC::toString(u_mean) << std::endl;
    //    std::cout << DSC::toString(u_i) << "   " << DSC::toString(u_j) << std::endl;

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(jacobian_eigen, jacobian_neg, jacobian_pos);
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector< EigenMatrixType > jacobian,
                           FluxJacobianRangeType& jacobian_neg,
                           FluxJacobianRangeType& jacobian_pos) const
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
        assert(std::abs(eigenvalues[ii].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      }
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real()*diag_jacobian_neg_tmp.backend()*eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real()*diag_jacobian_pos_tmp.backend()*eigenvectors_inverse.real());
      jacobian_neg[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_neg_eigen));
      jacobian_pos[ii] = DSC::fromString< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange > >(DSC::toString(jacobian_pos_eigen));
    }
#else
    static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const BoundaryValueFunctionType& boundary_values_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class Dirichlet

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
typename internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >::FluxJacobianRangeType Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp >::jacobian_neg_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
typename internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >::FluxJacobianRangeType Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp >::jacobian_pos_(0);

template < class LocalizableFunctionImp, class BoundaryValueFunctionImp >
bool Dirichlet< LocalizableFunctionImp, BoundaryValueFunctionImp >::jacobians_constructed_(false);


} // namespace Godunov
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_GODUNOV_HH
