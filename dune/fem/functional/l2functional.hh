
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

    L2Functional( const InducingFunctionType& inducingFunction )
      : inducingFunction_( inducingFunction )
    {
    }

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
      typedef typename EntityType::Geometry
        EntityGeometryType;
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
        const EntityGeometryType& entityGeometry = entity.geometry();

        // quadrature
        typedef CachingQuadrature< GridPartType, 0 >
          EntityQuadratureType;
        typedef typename EntityQuadratureType::CoordinateType
          EntityCoordinateType;
        const int quadratureOrder = (2 * discreteFunctionSpace.order() + 1);
        const EntityQuadratureType entityQuadrature( entity, quadratureOrder );
        const unsigned int numberOfQuadraturePoints = entityQuadrature.nop();

        // local function and basefunction set
        const LocalFunctionType& localFunction = discreteFunction.localFunction( entity );
        const BaseFunctionSetType baseFunctionSet = discreteFunctionSpace.baseFunctionSet( entity );
        const unsigned numberOfLocalDoFs = baseFunctionSet.numBaseFunctions();

        // local DoF and functional vector
        typedef LocalDoFVector< RangeFieldType >
          LocalDoFVectorType;
        const LocalDoFVectorType localDoFVector( localFunction );
        LocalDoFVectorType localFunctionalVector( numberOfLocalDoFs );

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
          localFunctionalVector[localDoF] = localFunctionalValue;

        } // done loop over all local DoFs

        // compute product
        ret += localDoFVector * localFunctionalVector;

      } // done gridwalk

      return ret;
    }

  private:

    const InducingFunctionType& inducingFunction_;

  }; // end class L2Functional

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
