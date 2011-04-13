/**
  * \file   integrator.hh
  * \brief  Contains various integrator classes.
  **/

#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/localoperation/interface.hh>

namespace Dune
{

namespace Functionals
{

namespace LocalOperation
{

namespace Integrator
{

/**
  * \brief  This class carries out a codim 0 integration.
  *
  *         Given a local operation \f$F\f$, which operates on one local function \f$u\f$ or on two local functions
  *         \f$u\f$, \f$v\f$, this class computes the integral
  *         \f$\int_{E} F(u) \text{dx}\f$ or \f$\int_{E} F(u,v) \text{dx}\f$ by a quadrature rule, where \f$E\f$ is the
  *         codim 0 entity, on which the local functions are defined.
  *
  * \tparam FunctionSpaceImp
  *         Type of the function space, the local functions live in, i.e. Dune::FunctionSpaceInterface.
  *
  * \tparam LocalOperationImp
  *         Type of the underlying local operation, i.e. Dune::Functionals::LocalOperation::Interface.
  **/
template< class FunctionSpaceImp, class LocalOperationImp >
class Codim0
  : public Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
{
public:

  typedef typename Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
    BaseType;

  typedef LocalOperationImp
    LocalOperationType;

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

private:

  // This Functor represents F(u) and is needed for the integrate() method.
  template< class LocalOperationType,  class LocalTestFunctionType >
  class SingleFunctor
  {
  public:

    typedef typename LocalTestFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename LocalTestFunctionType::EntityType
      EntityType;
    typedef typename LocalTestFunctionType::RangeFieldType
      RangeFieldType;

    SingleFunctor(  const LocalOperationType& localOperation,
                    const LocalTestFunctionType& localTestFunction )
      : localOperation_( localOperation ),
        localTestFunction_( localTestFunction )
    {
    }

    const EntityType& entity() const
    {
      return localTestFunction_.entity();
    }

    const unsigned int order() const
    {
      return localTestFunction_.order();
    }

    template< class LocalPointType >
    const RangeFieldType evaluate( const LocalPointType& x ) const
    {
      return localOperation_.evaluate( localTestFunction_, x );
    }

  private:

    const LocalOperationType& localOperation_;
    const LocalTestFunctionType& localTestFunction_;

  };

  // This Functor represents F(u,v) and is needed for the integrate() method.
  template< class LocalOperationType, class LocalAnsatzFunctionType, class LocalTestFunctionType >
  class DoubleFunctor
  {
  public:

    typedef typename LocalAnsatzFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename LocalAnsatzFunctionType::EntityType
      EntityType;
    typedef typename LocalAnsatzFunctionType::RangeFieldType
      RangeFieldType;

    DoubleFunctor(  const LocalOperationType& localOperation,
                    const LocalAnsatzFunctionType& localAnsatzFunction,
                    const LocalTestFunctionType& localTestFunction )
      : localOperation_( localOperation ),
        localAnsatzFunction_( localAnsatzFunction ),
        localTestFunction_( localTestFunction )
    {
    }

    const EntityType& entity() const
    {
      return localAnsatzFunction_.entity();
    }

    const unsigned int order() const
    {
      const unsigned int ansatzOrder = localAnsatzFunction_.order();
      const unsigned int testOrder = localTestFunction_.order();

      return std::max( ansatzOrder, testOrder );
    }

    template< class LocalPointType >
    const RangeFieldType evaluate( const LocalPointType& x ) const
    {
      return localOperation_.evaluate( localAnsatzFunction_, localTestFunction_, x );
    }

  private:

    const LocalOperationType& localOperation_;
    const LocalAnsatzFunctionType& localAnsatzFunction_;
    const LocalTestFunctionType& localTestFunction_;

  };

public:

  /**
    * \brief      Constructor storing the local operation.
    * \param[in]  localOperation
    *             The local Operation \f$F\f$.
    **/
  Codim0( const LocalOperationType& localOperation )
    : BaseType(),
      localOperation_( localOperation )
  {
    std::cout << "VolumeIntegrator::VolumeIntegrator()" << std::endl;
  }

  /**
    * \brief  Returns the local operation.
    * \return \f$F\f$
    **/
  LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  /**
    * \brief      Computes the integral \f$\int_{E} F(u) \text{dx}\f$.
    * \tparam     LocalFunctionType
    *             Type of the local function \f$u\f$, i.e. Dune::LocalFunction.
    * \param[in]  localTestFunction
    *             The local function \f$u\f$.
    * \return     The integral \f$\int_{E} F(u) \text{dx}\f$.
    **/
  template< class LocalTestFunctionType >
  RangeFieldType operate( const LocalTestFunctionType& localTestFunction ) const
  {
    typedef SingleFunctor< LocalOperationType, LocalTestFunctionType >
      FunctorType;

    return integrate( FunctorType( localOperation_, localTestFunction ) );
  }

  /**
    * \brief      Computes the integral \f$\int_{E} F(u,v) \text{dx}\f$.
    * \tparam     LocalAnsatzFunctionType
    *             Type of the local function \f$u\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalTestFunctionType
    *             Type of the local function \f$v\f$, i.e. Dune::LocalFunction.
    * \param[in]  localAnsatzFunction
    *             The local function \f$u\f$.
    * \param[in]  localTestFunction
    *             The local function \f$v\f$.
    * \return     The integral \f$\int_{E} F(u,v) \text{dx}\f$.
    **/
  template< class LocalAnsatzFunctionType, class LocalTestFunctionType >
  RangeFieldType operate( const LocalAnsatzFunctionType& localAnsatzFunction,
                          const LocalTestFunctionType& localTestFunction ) const
  {
    typedef DoubleFunctor< LocalOperationType, LocalAnsatzFunctionType, LocalTestFunctionType >
      FunctorType;

    return integrate( FunctorType( localOperation_, localAnsatzFunction, localTestFunction ) );
  }

private:

  // The actual integration is done here. The argument is either a SingleFunctor or a DoubleFunctor, so we can target
  // the integration of F(u) and F(u,v) with a single method.
  template< class FunctorType >
  const RangeFieldType integrate( const FunctorType& functor ) const
  {
    // init return value
    RangeFieldType ret = 0.0;

    // some types we will need
    typedef typename FunctorType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType
      GridPartType;

    typedef typename FunctorType::EntityType
      EntityType;

    typedef typename EntityType::Geometry
      EntityGeometryType;

    typedef CachingQuadrature< GridPartType, 0 >
      VolumeQuadratureType;

    // entity and geometry
    const EntityType& entity = functor.entity();
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    const unsigned int quadratureOrder = functor.order();
    const VolumeQuadratureType volumeQuadrature( entity, quadratureOrder );
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    for( unsigned int q = 0; q < numberOfQuadraturePoints; ++q )
    {
      // local coordinate
      const DomainType x = volumeQuadrature.point( q );

      // integration factors
      const double integrationFactor = entityGeometry.integrationElement( x );
      const double quadratureWeight = volumeQuadrature.weight( q );

      // evaluate the local operation
      const RangeFieldType localOperationEvalauted = functor.evaluate( x );

      // compute integral
      ret += integrationFactor * quadratureWeight * localOperationEvalauted;
    }

    // return
    return ret;
  }

  const LocalOperationType localOperation_;

}; // end class Codim0

} // end namespace Integrator

} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH
