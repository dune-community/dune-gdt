
#ifndef DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH
#define DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH

// dune common includes
#include <dune/common/fvector.hh>

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune fem-functionals includes
#include <dune/fem/common/localbasefunction.hh>
#include <dune/fem/common/localvector.hh>
#include <dune/fem/localoperation/interface.hh>

namespace Dune {

namespace Functionals {


namespace Functional {

/**
  * \brief      This class represents a finite element functional.
  *
  * \attention  This class is under construction!
  *
  * \todo       Doc me, please!
  **/
template <class DiscreteFunctionSpaceImp, class LocalOperationImp>
class FiniteElement
{
public:
  //! Type of the discrete function space.
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  //! Type of the local operation type.
  typedef LocalOperationImp LocalOperationType;

  //! Type of the analytical function space Type.
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! Type of domain vector (using type of domain field) has a Dune::FieldVector type interface.
  typedef typename FunctionSpaceType::DomainType DomainType;

  //! Intrinsic type used for values in the domain field (usually a double).
  typedef typename FunctionSpaceType::RangeFieldType DomainFieldType;

  //! Type of range vector (using type of range field) has a Dune::FieldVector type interface.
  typedef typename FunctionSpaceType::RangeType RangeType;

  //! Intrinsic type used for values in the range field (usually a double).
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  //! Type of the local vector.
  typedef Dune::Functionals::Common::LocalVector<RangeFieldType> LocalVectorType;

  //! Type of the basis function set.
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  //! Type of the grid part.
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  //! Type of the iterator iterating over entities.
  typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;

  //! Type of entities.
  typedef typename EntityIteratorType::Entity EntityType;

  //! Type of geometry of entities.
  typedef typename EntityType::Geometry EntityGeometryType;

  //! Type of quadrature on the grid part.
  typedef CachingQuadrature<GridPartType, 0> EntityQuadratureType;

  //! Type of the local basis function provider.
  typedef Dune::Functionals::Common::LocalBaseFunctionProvider<DiscreteFunctionSpaceType> LocalBaseFunctionProviderType;

  //! Type of local basis functions.
  typedef typename LocalBaseFunctionProviderType::LocalBaseFunctionType LocalBaseFunctionType;

  /**
   * @brief Constructor.
   *
   * @param discreteFunctionSpace The discrete function space.
   * @param localOperation The local operation.
   */
  FiniteElement(const DiscreteFunctionSpaceType& discreteFunctionSpace, const LocalOperationType& localOperation)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , localOperation_(localOperation)
    , localBaseFunctionProvider_(discreteFunctionSpace)
  {
  }

  /**
    * \brief  Copy constructor
    *
    * @param other Another finite element.
    **/
  FiniteElement(const FiniteElement& other)
    : discreteFunctionSpace_(other.space())
    , localOperation_(other.localOperation())
    , localBaseFunctionProvider_(discreteFunctionSpace_)
  {
  }

  /**
    * \brief  Destructor
    *         Does nothing.
    **/
  ~FiniteElement()
  {
  }

  /**
   * @brief Returns a reference to the discrete function space.
   *
   * @return A Reference to the discrete function space.
   */
  const DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

  /**
   * @brief Returns a local operation.
   *
   * @return A local operation.
   */
  const LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  /**
   *  @todo Doc me, please!
   */
  const LocalVectorType applyLocal(const EntityType& entity) const
  {
    // basefunctionset
    const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace_.baseFunctionSet(entity);
    const unsigned numberOfLocalDoFs          = baseFunctionSet.numBaseFunctions();

    // init return vector
    LocalVectorType ret(numberOfLocalDoFs);

    // geometry
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    const unsigned quadratureOrder = discreteFunctionSpace_.order();
    const EntityQuadratureType entityQuadrature(entity, quadratureOrder);
    const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

    // do loop over all local DoFs
    for (unsigned int i = 0; i < numberOfLocalDoFs; ++i) {
      // value of the functional, applied to the local basefunction
      RangeFieldType localFunctionalValue = 0.0;

      // get local basefunctions
      const LocalBaseFunctionType phi_i = localBaseFunctionProvider_.provide(entity, i);

      // do walk over quadrature points
      for (unsigned int quadraturePoint = 0; quadraturePoint < numberOfQuadraturePoints; ++quadraturePoint) {
        // coordinates
        const DomainType x = entityQuadrature.point(quadraturePoint);

        // integration factors
        const double integrationFactor = entityGeometry.integrationElement(x);
        const double quadratureWeight  = entityQuadrature.weight(quadraturePoint);

        // apply local operation
        const double localOperationEvaluated = localOperation_.operate(phi_i, x);

        // compute integral
        localFunctionalValue += integrationFactor * quadratureWeight * localOperationEvaluated;

      } // done walk over quadrature points

      // set local vector
      ret[i] = localFunctionalValue;

    } // done loop over all local DoFs

    return ret;

  } // end applyLocal()

private:
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const LocalOperationType localOperation_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class FiniteElement

/**
 * @brief This class represents ???
 *
 * @todo Doc me, please!
 *
 * @tparam DiscreteFunctionSpaceImp
 * @tparam LocalOperationImp
 */
template <class DiscreteFunctionSpaceImp, class LocalOperationImp>
class FiniteElementLOP
{
public:
  //! Type of the discrete function space.
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  //! Type of the analytical function space.
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  //! Type of the local operation.
  typedef LocalOperationImp LocalOperationType;

  //! Type ot the entity.
  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  //! Intrinsic type used for values in the range field (usually a double).
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  //! Type of the local dof vector.
  typedef Dune::Functionals::Common::LocalVector<RangeFieldType> LocalVectorType;

  //! Type of the local basis function provider.
  typedef Dune::Functionals::Common::LocalBaseFunctionProvider<DiscreteFunctionSpaceType> LocalBaseFunctionProviderType;

  //! Type of the local basis functions.
  typedef typename LocalBaseFunctionProviderType::LocalBaseFunctionType LocalBaseFunctionType;

  /**
   * @brief Constructor.
   *
   * @param discreteFunctionSpace The discrete function space.
   * @param localOperation The local operation.
   */
  FiniteElementLOP(const DiscreteFunctionSpaceType& discreteFunctionSpace, const LocalOperationType& localOperation)
    : discreteFunctionSpace_(discreteFunctionSpace)
    , localOperation_(localOperation)
    , localBaseFunctionProvider_(discreteFunctionSpace)
  {
  }

  /**
    * \brief  Destructor
    *         Does nothing.
    **/
  ~FiniteElementLOP()
  {
  }

  /**
   * @brief Returns a reference to the discrete function space.
   *
   * @return A Reference to the discrete function space.
   */
  const DiscreteFunctionSpaceType& space() const
  {
    return discreteFunctionSpace_;
  }

  /**
   * @brief Returns a local operation.
   *
   * @return A local operation.
   */
  LocalOperationType localOperation() const
  {
    return localOperation_;
  }

  /**
   * @todo Doc me, please!
   */
  LocalVectorType applyLocal(const EntityType& entity) const
  {
    const unsigned int numberOfLocalDoFs = discreteFunctionSpace_.baseFunctionSet(entity).numBaseFunctions();

    LocalVectorType ret(numberOfLocalDoFs);

    // do loop over all local DoFs
    for (unsigned int i = 0; i < numberOfLocalDoFs; ++i) {
      const LocalBaseFunctionType localBaseFunction_i = localBaseFunctionProvider_.provide(entity, i);

      ret[i] = localOperation_.operate(localBaseFunction_i);

    } // done loop over all local DoFs

    return ret;

  } // end applyLocal()

private:
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const LocalOperationType localOperation_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class FiniteElementLOP

} // end namespace Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_FUNCTIONAL_FINITEELEMENT_HH
