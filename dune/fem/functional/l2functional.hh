
#ifndef DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
#define DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH

// dune common includes
#include <dune/common/fvector.hh>

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune fem-tools includes
#include "../../../tools/function/functiontools.hh" // should be removed in the end!

namespace Dune
{

namespace Functionals
{

  /**
    * \brief      This class represents a local DoF vector.
    *
    *             It is based upon std::vector and should be replaced by something clever in the future!
    *
    * \todo       Doc me, please!
    **/
  template< class ElementType >
  class LocalDoFVector
  {
  public:
    /**
      * \brief    Initializes an empty vector, according to the given size.
      **/
    LocalDoFVector( const unsigned size )
      :size_( size )
    {
      // resize
      storage_.resize( size );
    }

    /**
      * \brief    Initializes a DoF vector and sets its entries to the
      *           corresponding entries of the given localFunction.
      **/
    template< class LocalFunctionType >
    LocalDoFVector( const LocalFunctionType& localFunction )
      :size_( localFunction.numDofs() )
    {
      // resize
      storage_.resize( localFunction.numDofs() );

      // copy entries
      for(  unsigned ii = 0;
            ii < localFunction.numDofs();
            ++ii )
      {
        storage_[ii] = localFunction[ii];
      }
    }

    /**
      * \brief    Returns the size.
      */
    unsigned size() const
    {
      return size_;
    }

    /**
      * \brief    Random read and write access.
      **/
    ElementType& operator[]( const unsigned ii )
    {
      return storage_[ii];
    }

    /**
      * \brief    Random read access.
      **/
    const ElementType operator[]( const unsigned ii ) const
    {
      return storage_[ii];
    }

    /**
      * \brief    Scalar product of two local DoF vectors of same type.
      **/
    ElementType operator*( const LocalDoFVector< ElementType >& other ) const
    {
      assert( size_ == other.size() );
      ElementType result = 0.0;

      for(  unsigned ii = 0;
            ii < size_;
            ++ii )
      {
        result += storage_[ii] * other[ii];
      }

      return result;
    }

  private:
    std::vector< ElementType > storage_;
    const unsigned size_;

  };

  /**
    * \brief      This class represents an L2 functional.
    *
    * \attention  This class is under construction!
    *
    * \todo       Doc me, please!
    **/
  template< class InducingFunctionImp >
  class L2Functional
  {
  public:

    typedef InducingFunctionImp
      InducingFunctionType;

    typedef typename InducingFunctionType::RangeFieldType
      RangeFieldType;

    typedef LocalDoFVector< RangeFieldType >
      LocalDoFVectorType;

    L2Functional( const InducingFunctionType& inducingFunction )
      : inducingFunction_( inducingFunction )
    {
    }

    /**
      * \brief      This function represents the application of the functional to a
      *             discrete function.
      *
      * \todo       Doc me, please!
      **/
    template< class DiscreteFunctionType >
    RangeFieldType operator()( const DiscreteFunctionType& discreteFunction ) const
    {
      RangeFieldType ret = 0.0;

      // some types we will need
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
        BaseFunctionSetType;
      typedef typename DiscreteFunctionSpaceType::GridPartType
        GridPartType;
      typedef typename DiscreteFunctionSpaceType::IteratorType
        EntityIteratorType;
      typedef typename EntityIteratorType::Entity
        EntityType;
      typedef typename DiscreteFunctionType::RangeType
        RangeType;
      typedef typename DiscreteFunctionType::LocalFunctionType
        LocalFunctionType;

      // some things we will need
      const DiscreteFunctionSpaceType& discreteFunctionSpace = discreteFunction.space();

      // do gridwalk
      const EntityIteratorType BehindLastEntityIterator = discreteFunctionSpace.end();
      for ( EntityIteratorType entityIterator = discreteFunctionSpace.begin();
            entityIterator != BehindLastEntityIterator;
            ++entityIterator )
      {
        // entity and geometry
        const EntityType& entity = *entityIterator;

        // local function and basefunction set
        const LocalFunctionType& localFunction = discreteFunction.localFunction( entity );
        const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace.baseFunctionSet( entity );

        // local DoF and functional vector
        const LocalDoFVectorType localDoFVector( localFunction );
        const LocalDoFVectorType localFunctionalVector = applyLocal( discreteFunctionSpace, entity, baseFunctionSet );

        // compute product
        ret += localDoFVector * localFunctionalVector;

      } // done gridwalk

      return ret;

    } // end operator()

    /**
      * \brief      This function respresents the local application of the functional
      *             to a local basefunctionset.
      *
      * \todo       Doc me, please!
      **/
    template< class DiscreteFunctionSpaceType, class EntityType, class BaseFunctionSetType >
    LocalDoFVectorType applyLocal( const DiscreteFunctionSpaceType& discreteFunctionSpace, const EntityType& entity, const BaseFunctionSetType& baseFunctionSet ) const
    {
      // some types we will need
      typedef typename DiscreteFunctionSpaceType::GridPartType
        GridPartType;
      typedef typename DiscreteFunctionSpaceType::IteratorType
        EntityIteratorType;
      typedef typename EntityIteratorType::Entity
        EntityType;
      typedef typename EntityType::Geometry
        EntityGeometryType;
      typedef typename BaseFunctionSetType::RangeType
        RangeType;

      // init return vector
      const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();
      LocalDoFVectorType ret( numberOfLocalDoFs );

      // geometry
      const EntityGeometryType& entityGeometry = entity.geometry();

      // quadrature
      typedef CachingQuadrature< GridPartType, 0 >
        EntityQuadratureType;
      typedef typename EntityQuadratureType::CoordinateType
        EntityCoordinateType;
      const unsigned quadratureOrder = (2 * discreteFunctionSpace.order() + 1);
      const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
      const unsigned numberOfQuadraturePoints = entityQuadrature.nop();

      // do loop over all local DoFs
      for(  unsigned int localDoF = 0;
            localDoF < numberOfLocalDoFs;
            ++localDoF )
      {
        // value of the L2 functional, applied to the local basefunction, associated with the local DoF
        RangeFieldType localFunctionalValue = 0.0;
        // do walk over quadrature points
        for(  unsigned int quadraturePoint = 0;
              quadraturePoint < numberOfQuadraturePoints;
              ++quadraturePoint )
        {
          // coordinates
          const EntityCoordinateType xReferenceElement = entityQuadrature.point( quadraturePoint );
          const EntityCoordinateType xWorld = entityGeometry.global( xReferenceElement );

          // integration factors
          const double integrationFactor = entityGeometry.integrationElement( xReferenceElement );
          const double quadratureWeight = entityQuadrature.weight( quadraturePoint );

          // evaluate function and basefunction
          RangeType functionValue = 0.0;
          inducingFunction_.evaluate( xWorld, functionValue );
          RangeType baseFunctionValue = 0.0;
          baseFunctionSet.evaluate( localDoF, xReferenceElement, baseFunctionValue );

          // compute integral
          localFunctionalValue += integrationFactor * quadratureWeight * functionValue * baseFunctionValue;

        } // done walk over quadrature points

        // set local vector
        ret[localDoF] = localFunctionalValue;

      } // done loop over all local DoFs

      return ret;

    } // end applyLocal()

  private:

    const InducingFunctionType& inducingFunction_;

  }; // end class L2Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
