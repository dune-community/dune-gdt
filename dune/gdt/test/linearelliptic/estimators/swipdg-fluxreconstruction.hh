// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_ESTIMATORS_SWIPDG_HH
#define DUNE_HDD_LINEARELLIPTIC_ESTIMATORS_SWIPDG_HH

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/grid/walker/functors.hh>
#include <dune/stuff/playground/functions/ESV2007.hh>

#include <dune/pymor/common/exceptions.hh>
#include <dune/pymor/parameters/base.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/playground/localevaluation/ESV2007.hh>
#include <dune/gdt/playground/operators/fluxreconstruction.hh>
#include <dune/gdt/playground/spaces/finitevolume/default.hh>

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Estimators {
namespace internal {
namespace SWIPDG {


static const size_t over_integrate = 2;


class LocalNonconformityESV2007Base
{
public:
  static std::string id() { return "eta_NC_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class LocalNonconformityESV2007
  : public LocalNonconformityESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local nonconformity estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType >
class LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalNonconformityESV2007Base
  , public Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
    ThisType;
  typedef Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType   RangeFieldType;

private:
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
  typedef typename ConstDiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Elliptic< DiffusionFactorType,
                                                                              DiffusionTensorType > > LocalOperatorType;
  typedef Stuff::Common::TmpMatricesStorage< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem, const Pymor::Parameter& mu_bar)
  {
    if (mu_bar.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu_bar is of type " << mu_bar.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric diffusion_tensor!");
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_bar") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_bar'!");
    const Pymor::Parameter mu_bar = problem.parametric() ? parameters.at("mu_bar") : Pymor::Parameter();
    ThisType estimator(space, vector, problem, mu_bar);
    Stuff::Grid::Walker< GridViewType > grid_walker(space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalNonconformityESV2007(const SpaceType& space,
                            const VectorType& vector,
                            const ProblemType& problem,
                            const Pymor::Parameter mu_bar = Pymor::Parameter())
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem, mu_bar))
    , problem_mu_bar_(problem_.with_mu(mu_bar))
    , discrete_solution_(space_, vector_)
    , oswald_interpolation_(space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(discrete_solution_ - oswald_interpolation_))
    , local_operator_(over_integrate,
                      *problem_mu_bar_->diffusion_factor()->affine_part(),
                      *problem_.diffusion_tensor()->affine_part())
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , prepared_(false)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    if (!prepared_) {
      const GDT::Operators::OswaldInterpolation< GridViewType > oswald_interpolation_operator(space_.grid_view());
      oswald_interpolation_operator.apply(discrete_solution_, oswald_interpolation_);
      result_ = 0.0;
      prepared_ = true;
    }
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_difference = difference_->local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const SpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_bar_;
  const ConstDiscreteFunctionType discrete_solution_;
  DiscreteFunctionType oswald_interpolation_;
  std::unique_ptr< const DifferenceType > difference_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
  bool prepared_;
public:
  RangeFieldType result_;
}; // class LocalNonconformityESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class LocalResidualESV2007Base
{
public:
  static std::string id() { return "eta_R_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class LocalResidualESV2007
  : public LocalResidualESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local residual estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType >
class LocalResidualESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalResidualESV2007Base
  , public Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalResidualESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
      ThisType;
  typedef Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

private:
  typedef GDT::Spaces::FiniteVolume::Default< GridViewType, RangeFieldType, 1, 1 > P0SpaceType;
  typedef GDT::DiscreteFunction< P0SpaceType, VectorType > DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef typename Stuff::Functions::ESV2007::Cutoff< DiffusionFactorType, DiffusionTensorType > CutoffFunctionType;
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Product< CutoffFunctionType > > LocalOperatorType;
  typedef DSC::TmpMatricesStorage< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem)
  {
    if (problem.parametric())
      DUNE_THROW(NotImplemented, "Not implemented yet for parametric problems!");
    assert(problem.diffusion_factor()->has_affine_part());
    assert(problem.diffusion_tensor()->has_affine_part());
    assert(problem.force()->has_affine_part());
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space, const VectorType& /*vector*/, const ProblemType& problem)
  {
    ThisType estimator(space, problem);
    Stuff::Grid::Walker< GridViewType > grid_walker(space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalResidualESV2007(const SpaceType& space, const ProblemType& problem)
    : space_(space)
    , problem_(assert_problem(problem))
    , p0_space_(space_.grid_view())
    , p0_force_(p0_space_)
    , difference_(Stuff::Common::make_unique< DifferenceType >(*problem_.force()->affine_part() - p0_force_))
    , cutoff_function_(*problem_.diffusion_factor()->affine_part(),
                       *problem_.diffusion_tensor()->affine_part())
    , local_operator_(over_integrate, cutoff_function_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , prepared_(false)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    if (!prepared_) {
      const GDT::Operators::Projection< GridViewType > projection_operator(space_.grid_view(), over_integrate);
      projection_operator.apply(*problem_.force()->affine_part(), p0_force_);
      result_ = 0.0;
      prepared_ = true;
    }
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_difference = difference_->local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const SpaceType& space_;
  const ProblemType& problem_;
  const P0SpaceType p0_space_;
  DiscreteFunctionType p0_force_;
  std::unique_ptr< const DifferenceType > difference_;
  const CutoffFunctionType cutoff_function_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
  bool prepared_;
public:
  RangeFieldType result_;
}; // class LocalResidualESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class LocalResidualESV2007StarBase
{
public:
  static std::string id() { return "eta_R_ESV2007_*"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class LocalResidualESV2007Star
  : public LocalResidualESV2007StarBase
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local residual estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType >
class LocalResidualESV2007Star< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalResidualESV2007StarBase
  , public Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalResidualESV2007Star< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
      ThisType;
  typedef Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static const unsigned int dimDomain = SpaceType::dimDomain;

private:
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::Spaces::RaviartThomas::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
  typedef GDT::DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;
  typedef typename RTN0DiscreteFunctionType::DivergenceType DivergenceType;
  typedef typename DivergenceType::DifferenceType DifferenceType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef typename Stuff::Functions::ESV2007::Cutoff< DiffusionFactorType, DiffusionTensorType > CutoffFunctionType;
  typedef GDT::LocalOperator::Codim0Integral< GDT::LocalEvaluation::Product< CutoffFunctionType > > LocalOperatorType;
  typedef DSC::TmpMatricesStorage< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem, const Pymor::Parameter& mu)
  {
    if (mu.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu is of type " << mu.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric diffusion_tensor!");
    if (problem.force()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric force!");
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (problem.diffusion_factor()->parametric() && parameters.find("mu") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu'!");
    const Pymor::Parameter mu = problem.parametric() ? parameters.at("mu") : Pymor::Parameter();
    ThisType estimator(space, vector, problem, mu);
    Stuff::Grid::Walker< GridViewType > grid_walker(space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalResidualESV2007Star(const SpaceType& space,
                           const VectorType& vector,
                           const ProblemType& problem,
                           const Pymor::Parameter mu = Pymor::Parameter())
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem, mu))
    , problem_mu_(problem_.with_mu(mu))
    , discrete_solution_(space_, vector_)
    , rtn0_space_(space_.grid_view())
    , diffusive_flux_(rtn0_space_)
    , divergence_(diffusive_flux_.divergence())
    , difference_(*problem_.force()->affine_part() - divergence_)
    , cutoff_function_(*problem_.diffusion_factor()->affine_part(),
                       *problem_.diffusion_tensor()->affine_part())
    , local_operator_(over_integrate, cutoff_function_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , prepared_(false)
    , result_(0.0)
  {}

  virtual ~LocalResidualESV2007Star() = default;

  virtual void prepare()
  {
    if (!prepared_) {
      const GDT::Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionFactorType, DiffusionTensorType >
        diffusive_flux_reconstruction(space_.grid_view(),
                                      *problem_mu_->diffusion_factor()->affine_part(),
                                      *problem_.diffusion_tensor()->affine_part(),
                                      over_integrate);
      diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
      result_ = 0.0;
      prepared_ = true;
    }
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_difference = difference_.local_function(entity);
    local_operator_.apply(*local_difference,
                          *local_difference,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const SpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const std::shared_ptr< typename ProblemType::NonparametricType > problem_mu_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const DivergenceType divergence_;
  const DifferenceType difference_;
  const CutoffFunctionType cutoff_function_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
  bool prepared_;
public:
  RangeFieldType result_;
}; // class LocalResidualESV2007Star< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class LocalDiffusiveFluxESV2007Base
{
public:
  static std::string id() { return "eta_DF_ESV2007"; }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class LocalDiffusiveFluxESV2007
  : public LocalDiffusiveFluxESV2007Base
{
public:
  static const bool available = false;
};

#if HAVE_ALUGRID

/**
 *  \brief computes the local diffusive flux estimator as defined in ESV2007
 */
template< class SpaceType, class VectorType, class ProblemType >
class LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public LocalDiffusiveFluxESV2007Base
  , public Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType >
{
  typedef LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
      ThisType;
  typedef Stuff::Grid::Functor::Codim0< typename SpaceType::GridViewType > FunctorBaseType;
public:
  static const bool available = true;

  typedef std::map< std::string, Pymor::Parameter > ParametersMapType;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename FunctorBaseType::EntityType   EntityType;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static const unsigned int dimDomain = SpaceType::dimDomain;

private:
  typedef GDT::ConstDiscreteFunction< SpaceType, VectorType > ConstDiscreteFunctionType;
  typedef GDT::Spaces::RaviartThomas::PdelabBased< GridViewType, 0, RangeFieldType, dimDomain > RTN0SpaceType;
  typedef GDT::DiscreteFunction< RTN0SpaceType, VectorType > RTN0DiscreteFunctionType;

  typedef typename ProblemType::DiffusionFactorType::NonparametricType DiffusionFactorType;
  typedef typename ProblemType::DiffusionTensorType::NonparametricType DiffusionTensorType;

  typedef GDT::LocalOperator::Codim0Integral<
      GDT::LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionFactorType,
                                                            RTN0DiscreteFunctionType,
                                                            DiffusionTensorType > > LocalOperatorType;
  typedef DSC::TmpMatricesStorage< RangeFieldType > TmpStorageProviderType;

  static const ProblemType& assert_problem(const ProblemType& problem,
                                           const Pymor::Parameter& mu,
                                           const Pymor::Parameter& mu_hat)
  {
    if (mu.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu is of type " << mu.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (mu_hat.type() != problem.parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                 "Given mu_hat is of type " << mu_hat.type() << " and should be of type " << problem.parameter_type()
                 << "!");
    if (problem.diffusion_tensor()->parametric())
      DUNE_THROW(NotImplemented, "Not implemented for parametric diffusion_tensor!");
    return problem;
  } // ... assert_problem(...)

public:
  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const ParametersMapType parameters = ParametersMapType())
  {
    if (problem.diffusion_factor()->parametric() && parameters.find("mu") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu'!");
    if (problem.diffusion_factor()->parametric() && parameters.find("mu_hat") == parameters.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, "Given parameters are missing 'mu_hat'!");
    const Pymor::Parameter mu =     problem.parametric() ? parameters.at("mu")     : Pymor::Parameter();
    const Pymor::Parameter mu_hat = problem.parametric() ? parameters.at("mu_hat") : Pymor::Parameter();
    ThisType estimator(space, vector, problem, mu, mu_hat);
    Stuff::Grid::Walker< GridViewType > grid_walker(space.grid_view());
    grid_walker.add(estimator);
    grid_walker.walk();
    return std::sqrt(estimator.result_);
  } // ... estimate(...)

  LocalDiffusiveFluxESV2007(const SpaceType& space,
                            const VectorType& vector,
                            const ProblemType& problem,
                            const Pymor::Parameter mu = Pymor::Parameter(),
                            const Pymor::Parameter mu_hat = Pymor::Parameter())
    : space_(space)
    , vector_(vector)
    , problem_(assert_problem(problem, mu, mu_hat))
    , problem_mu_(problem_.with_mu(mu))
    , problem_mu_hat_(problem.with_mu(mu_hat))
    , discrete_solution_(space_, vector_)
    , rtn0_space_(space.grid_view())
    , diffusive_flux_(rtn0_space_)
    , local_operator_(over_integrate,
                      *problem_mu_hat_->diffusion_factor()->affine_part(),
                      *problem_.diffusion_tensor()->affine_part(),
                      diffusive_flux_)
    , tmp_local_matrices_({1, local_operator_.numTmpObjectsRequired()}, 1, 1)
    , prepared_(false)
    , result_(0.0)
  {}

  virtual void prepare()
  {
    if (!prepared_) {
      const GDT::Operators::DiffusiveFluxReconstruction< GridViewType, DiffusionFactorType, DiffusionTensorType >
        diffusive_flux_reconstruction(space_.grid_view(),
                                      *problem_mu_->diffusion_factor()->affine_part(),
                                      *problem_.diffusion_tensor()->affine_part(),
                                      over_integrate);
      diffusive_flux_reconstruction.apply(discrete_solution_, diffusive_flux_);
      result_ = 0.0;
      prepared_ = true;
    }
  } // ... prepare(...)

  RangeFieldType compute_locally(const EntityType& entity)
  {
    const auto local_discrete_solution = discrete_solution_.local_function(entity);
    local_operator_.apply(*local_discrete_solution,
                          *local_discrete_solution,
                          tmp_local_matrices_.matrices()[0][0],
                          tmp_local_matrices_.matrices()[1]);
    assert(tmp_local_matrices_.matrices()[0][0].rows() >= 1);
    assert(tmp_local_matrices_.matrices()[0][0].cols() >= 1);
    return tmp_local_matrices_.matrices()[0][0][0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType &entity)
  {
    result_ += compute_locally(entity);
  }

private:
  const SpaceType& space_;
  const VectorType& vector_;
  const ProblemType& problem_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_;
  const std::shared_ptr< const typename ProblemType::NonparametricType > problem_mu_hat_;
  const ConstDiscreteFunctionType discrete_solution_;
  const RTN0SpaceType rtn0_space_;
  RTN0DiscreteFunctionType diffusive_flux_;
  const LocalOperatorType local_operator_;
  TmpStorageProviderType tmp_local_matrices_;
  bool prepared_;
public:
  RangeFieldType result_;
}; // class LocalDiffusiveFluxESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class ESV2007Base
{
public:
  static std::string id()
  {
    return "eta_ESV2007";
  }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class ESV2007
  : public ESV2007Base
{
public:
  static const bool available = false;
};


#if HAVE_ALUGRID

template< class SpaceType, class VectorType, class ProblemType >
class ESV2007< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public ESV2007Base
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;
public:
  static const bool available = true;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_squared(0.0);

    const auto& grid_view = space.grid_view();
    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_squared += eta_nc.compute_locally(entity)
                   + std::pow(std::sqrt(eta_r.compute_locally(entity)) + std::sqrt(eta_df.compute_locally(entity)), 2);
    }
    return std::sqrt(eta_squared);
  } // ... estimate(...)

  static Stuff::LA::CommonDenseVector< RangeFieldType > estimate_local(const SpaceType& space,
                                                                       const VectorType& vector,
                                                                       const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    const auto& grid_view = space.grid_view();
    Stuff::LA::CommonDenseVector< RangeFieldType >
        local_indicators(boost::numeric_cast< size_t >(grid_view.indexSet().size(0)), 0.0);
    RangeFieldType eta_squared = 0.0;

    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto index = grid_view.indexSet().index(entity);
      const RangeFieldType eta_t_squared
          = eta_nc.compute_locally(entity)
            + std::pow(std::sqrt(eta_r.compute_locally(entity)) + std::sqrt(eta_df.compute_locally(entity)), 2);
      local_indicators[index] = eta_t_squared;
      eta_squared += eta_t_squared;
    }
    for (auto& element : local_indicators)
      element /= eta_squared;
    return local_indicators;
  } // ... estimate_local(...)
}; // class ESV2007< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


class ESV2007AlternativeSummationBase
{
public:
  static std::string id()
  {
    return "eta_ESV2007_alt";
  }
};


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class ESV2007AlternativeSummation
  : public ESV2007AlternativeSummationBase
{
public:
  static const bool available = false;
};


#if HAVE_ALUGRID

template< class SpaceType, class VectorType, class ProblemType >
class ESV2007AlternativeSummation< SpaceType, VectorType, ProblemType, ALUGrid< 2, 2, simplex, conforming > >
  : public ESV2007AlternativeSummationBase
{
  typedef ALUGrid< 2, 2, simplex, conforming > GridType;
public:
  static const bool available = true;

  typedef typename ProblemType::RangeFieldType RangeFieldType;

  static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    RangeFieldType eta_nc_squared(0.0);
    RangeFieldType eta_r_squared(0.0);
    RangeFieldType eta_df_squared(0.0);

    const auto& grid_view = space.grid_view();
    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      eta_nc_squared += eta_nc.compute_locally(entity);
      eta_r_squared += eta_r.compute_locally(entity);
      eta_df_squared += eta_df.compute_locally(entity);
    }
    return std::sqrt(eta_nc_squared) + std::sqrt(eta_r_squared) + std::sqrt(eta_df_squared);
  } // ... estimate(...)

  static Stuff::LA::CommonDenseVector< RangeFieldType > estimate_local(const SpaceType& space,
                                                                       const VectorType& vector,
                                                                       const ProblemType& problem)
  {
    LocalNonconformityESV2007< SpaceType, VectorType, ProblemType, GridType > eta_nc(space, vector, problem);
    LocalResidualESV2007< SpaceType, VectorType, ProblemType, GridType >      eta_r(space, problem);
    LocalDiffusiveFluxESV2007< SpaceType, VectorType, ProblemType, GridType > eta_df(space, vector, problem);
    eta_nc.prepare();
    eta_r.prepare();
    eta_df.prepare();

    const auto grid_view = space.grid_view();
    Stuff::LA::CommonDenseVector< RangeFieldType >
        local_indicators(boost::numeric_cast< size_t >(grid_view.indexSet().size(0)), 0.0);
    RangeFieldType eta_nc_squared(0.0);
    RangeFieldType eta_r_squared(0.0);
    RangeFieldType eta_df_squared(0.0);

    const auto entity_it_end = grid_view.template end< 0 >();
    for (auto entity_it = grid_view.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto index = grid_view.indexSet().index(entity);
      const RangeFieldType eta_nc_t_squared = eta_nc.compute_locally(entity);
      const RangeFieldType eta_r_t_squared = eta_r.compute_locally(entity);
      const RangeFieldType eta_df_t_squared = eta_df.compute_locally(entity);
      eta_nc_squared += eta_nc_t_squared;
      eta_r_squared += eta_r_t_squared;
      eta_df_squared += eta_df_t_squared;
      local_indicators[index] = 3.0 *(eta_nc_t_squared + eta_r_t_squared + eta_df_t_squared);
    }
    const RangeFieldType eta_squared
        = std::pow(std::sqrt(eta_nc_squared) + std::sqrt(eta_r_squared) + std::sqrt(eta_df_squared), 2);
    for (auto& element : local_indicators)
      element /= eta_squared;
    return local_indicators;
  } // ... estimate_local(...)
}; // class ESV2007AlternativeSummation< ..., ALUGrid< 2, 2, simplex, conforming >, ... >

#endif // HAVE_ALUGRID


} // namespace SWIPDG
} // namespace internal


template< class SpaceType, class VectorType, class ProblemType, class GridType >
class SWIPDG
{
public:
  typedef typename ProblemType::RangeFieldType RangeFieldType;

private:
  template< class IndividualEstimator, bool available = false >
  class Caller
  {
  public:
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      return in;
    }

    static bool equals(const std::string& /*type*/)
    {
      return false;
    }

    static RangeFieldType estimate(const SpaceType& /*space*/,
                                   const VectorType& /*vector*/,
                                   const ProblemType& /*problem*/)
    {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
      return RangeFieldType(0);
    }

    static Stuff::LA::CommonDenseVector< RangeFieldType > estimate_local(const SpaceType& /*space*/,
                                                                         const VectorType& /*vector*/,
                                                                         const ProblemType& /*problem*/)
    {
      DUNE_THROW(Stuff::Exceptions::internal_error, "This should not happen!");
      return Stuff::LA::CommonDenseVector< RangeFieldType >();
    }
  }; // class Caller

  template< class IndividualEstimator >
  class Caller< IndividualEstimator, true >
  {
  public:
    static std::vector< std::string > append(std::vector< std::string > in)
    {
      in.push_back(IndividualEstimator::id());
      return in;
    }

    static bool equals(const std::string& type)
    {
      return IndividualEstimator::id() == type;
    }

    static RangeFieldType estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
    {
      return IndividualEstimator::estimate(space, vector, problem);
    }

    static Stuff::LA::CommonDenseVector< RangeFieldType > estimate_local(const SpaceType& space,
                                                                         const VectorType& vector,
                                                                         const ProblemType& problem)
    {
      return IndividualEstimator::estimate_local(space, vector, problem);
    }
  }; // class Caller< ..., true >

  template< class IndividualEstimator >
  static std::vector< std::string > call_append(std::vector< std::string > in)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::append(in);
  }

  template< class IndividualEstimator >
  static bool call_equals(const std::string& type)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::equals(type);
  }

  template< class IndividualEstimator >
  static RangeFieldType call_estimate(const SpaceType& space, const VectorType& vector, const ProblemType& problem)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::estimate(space, vector, problem);
  }

  template< class IndividualEstimator >
  static Stuff::LA::CommonDenseVector< RangeFieldType > call_estimate_local(const SpaceType& space,
                                                                            const VectorType& vector,
                                                                            const ProblemType& problem)
  {
    return Caller< IndividualEstimator, IndividualEstimator::available >::estimate_local(space, vector, problem);
  }

  typedef internal::SWIPDG::LocalNonconformityESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalNonconformityESV2007Type;
  typedef internal::SWIPDG::LocalResidualESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalResidualESV2007Type;
  typedef internal::SWIPDG::LocalResidualESV2007Star
      < SpaceType, VectorType, ProblemType, GridType >              LocalResidualESV2007StarType;
  typedef internal::SWIPDG::LocalDiffusiveFluxESV2007
      < SpaceType, VectorType, ProblemType, GridType >              LocalDiffusiveFluxESV2007Type;
  typedef internal::SWIPDG::ESV2007
      < SpaceType, VectorType, ProblemType, GridType >              ESV2007Type;
  typedef internal::SWIPDG::ESV2007AlternativeSummation
      < SpaceType, VectorType, ProblemType, GridType >              ESV2007AlternativeSummationType;

public:
  static std::vector< std::string > available()
  {
    std::vector< std::string > tmp;
    tmp = call_append< LocalNonconformityESV2007Type >(tmp);
    tmp = call_append< LocalResidualESV2007Type >(tmp);
    tmp = call_append< LocalResidualESV2007StarType >(tmp);
    tmp = call_append< LocalDiffusiveFluxESV2007Type >(tmp);
    tmp = call_append< ESV2007Type >(tmp);
    tmp = call_append< ESV2007AlternativeSummationType >(tmp);
    return tmp;
  } // ... available(...)

  static std::vector< std::string > available_local()
  {
    std::vector< std::string > tmp;
    tmp = call_append< ESV2007Type >(tmp);
    tmp = call_append< ESV2007AlternativeSummationType >(tmp);
    return tmp;
  } // ... available_local(...)

  static RangeFieldType estimate(const SpaceType& space,
                                 const VectorType& vector,
                                 const ProblemType& problem,
                                 const std::string type)
  {
    if (call_equals< LocalNonconformityESV2007Type >(type))
      return call_estimate< LocalNonconformityESV2007Type >(space, vector, problem);
    else if (call_equals< LocalResidualESV2007Type >(type))
      return call_estimate< LocalResidualESV2007Type >(space, vector, problem);
    else if (call_equals< LocalResidualESV2007StarType >(type))
      return call_estimate< LocalResidualESV2007StarType >(space, vector, problem);
    else if (call_equals< LocalDiffusiveFluxESV2007Type >(type))
      return call_estimate< LocalDiffusiveFluxESV2007Type >(space, vector, problem);
    else if (call_equals< ESV2007Type >(type))
      return call_estimate< ESV2007Type >(space, vector, problem);
    else if (call_equals< ESV2007AlternativeSummationType >(type))
      return call_estimate< ESV2007AlternativeSummationType >(space, vector, problem);
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Requested type '" << type << "' is not one of available()!");
  } // ... estimate(...)

  static Stuff::LA::CommonDenseVector< RangeFieldType > estimate_local(const SpaceType& space,
                                                                       const VectorType& vector,
                                                                       const ProblemType& problem,
                                                                       const std::string type)
  {
    if (call_equals< ESV2007Type >(type))
      return call_estimate_local< ESV2007Type >(space, vector, problem);
    else if (call_equals< ESV2007AlternativeSummationType >(type))
      return call_estimate_local< ESV2007AlternativeSummationType >(space, vector, problem);
    else
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Requested type '" << type << "' is not one of available_local()!");
  } // ... estimate_local(...)
}; // class SWIPDG


} // namespace Discretizations
} // namespace Estimators
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_ESTIMATORS_SWIPDG_HH
