
#ifndef DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
#define DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH

////- Dune includes
//#include <dune/common/fmatrix.hh>
//#include <dune/common/timer.hh>

//#include <dune/fem/storage/array.hh>
//#include <dune/fem/quadrature/quadrature.hh>
//#include <dune/fem/function/common/scalarproducts.hh>
//#include <dune/fem/operator/common/operator.hh>
//#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
//#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

//#include <dune/fem-howto/probleminterfaces.hh>

namespace Dune
{

  /**
    * \brief
    *
    * \todo doc me please
    **/
  template< class FunctionImp, class DiscreteFunctionImp >
  class LinearCodimZeroFunctional
  {
  public:

    typedef FunctionImp
      FunctionType;
    typedef DiscreteFunctionImp
      DiscreteFunctionType;

    typedef typename FunctionType::RangeType
      RangeType;
    typedef typename FunctionType::RangeFieldType
      RangeFieldType;


		typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType         DiscreteFunctionSpaceType;

    LinearCodimZeroFunctional( const FunctionType& function )
      : function_( function )
    {
    }

    RangeFieldType operator()( const DiscreteFunctionType &discreteFunction ) const
    {
      RangeFieldType ret( 0.0 );

      return ret;
    }

    template < class LocalFunctionType >
    void applyLocal( const LocalFunctionType &localFunction ) const
    {
    }

  private:
    const FunctionType &function_;

  };






  // L2 functional
  template< class Function, class DiscreteFunction >
  class IntegralFunctional
  {
  public:
    typedef Function FunctionType ;
    typedef typename FunctionType :: RangeType      RangeType ;
    typedef typename FunctionType :: RangeFieldType RangeFieldType ;

    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;

    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;

    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;

  protected:
    const FunctionType& function_;
    mutable std::vector< RangeType > functionValues_;

		const DiscreteFunctionSpaceType& dfSpace_;

  public:

		IntegralFunctional( const FunctionType &function, const DiscreteFunctionSpaceType& dfSpace )
      : function_( function ),
        functionValues_ (),
				dfSpace_(dfSpace)
    {
    }


    RangeFieldType operator () ( const DiscreteFunction& discreteFunction ) const
    {
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();

      Dune::TemporaryLocalFunction< DiscreteFunctionSpaceType >
        localFunctional( discreteFunctionSpace );

      const int quadratureOrder = (2 * discreteFunctionSpace.order() + 1);

      RangeFieldType functional = 0 ;
      // start grid traversal
      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin();
           it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        // add to global DoF vector
        LocalFunctionType localFunction = discreteFunction.localFunction( entity );

        // init local function
        localFunctional.init( entity );

        // apply function locally
        applyLocal( localFunctional, quadratureOrder );
        const int numDofs = localFunctional.numDofs();
        for(int l = 0; l < numDofs; ++l )
          functional += localFunctional[ l ] * localFunction[ l ];
      }

      return functional;
    }

    template <class LocalFunctionType>
    void applyLocal( LocalFunctionType& localFunctional,
                     const int quadOrder = -1 ) const
    {
      // clear values
      localFunctional.clear();

      const int quadratureOrder = ( quadOrder < 0 ) ?
        ( 2 * localFunctional.order() + 1 ) : quadOrder ;

      // We use a caching quadrature for codimension 0 entities
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

      // *it gives a reference to the current entity
      const EntityType &entity = localFunctional.entity();

      // obtain a reference to the entity's geometry
      const GeometryType &geometry = entity.geometry();

      QuadratureType quadrature( entity, quadratureOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      functionValues_.resize( numQuadraturePoints );
      for( size_t qP = 0; qP < numQuadraturePoints; ++qP )
      {
        // get integration element multiplied with quadrature weight
        const double factor =
          geometry.integrationElement( quadrature.point( qP ) ) *
          quadrature.weight( qP );

        // evaluate right hand side function
        function_.evaluate( geometry.global( quadrature.point( qP ) ), functionValues_[ qP ] );

        // apply factor
        functionValues_[ qP ] *= factor;
      }

      // add to local function (evaluates each base function)
      localFunctional.axpyQuadrature( quadrature, functionValues_ );
    }

    void algebraic( DiscreteFunction& discreteFunction, const int quadOrder = -1 ) const
    {
      const int quadratureOrder = ( quadOrder < 0 ) ? (2 * dfSpace_.order() + 1) : quadOrder;
      const DiscreteFunctionSpaceType &discreteFunctionSpace = dfSpace_;

			//DiscreteFunction discreteFunction("",dfSpace_);

      // set discreteFunction to zero
      discreteFunction.clear();

      Dune::TemporaryLocalFunction< DiscreteFunctionSpaceType >
        localFunctional( discreteFunctionSpace );

      // start grid traversal
      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin();
           it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        // init local function
        localFunctional.init( entity );

        // apply function locally
        this->applyLocal( localFunctional, quadratureOrder );

        // add to global DoF vector
        typedef typename DiscreteFunctionType :: LocalFunctionType  LocalFunctionType;
        LocalFunctionType dfLocalFunctional = discreteFunction.localFunction( entity );
        dfLocalFunctional += localFunctional;
      }

      // communicate data (for parallel runs)
      discreteFunction.communicate();
    }


  };

  // assembled functional
  template< class Functional >
  class AssembledFunctional
  {
  public:
    typedef Functional  FunctionalType ;
    typedef typename Functional :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;

    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;

    enum { dimension = GridType :: dimension };

    const DiscreteFunctionSpaceType& space_;
    const FunctionalType& functional_;
  public:
    AssembledFunctional(const DiscreteFunctionSpaceType& space,
                        const FunctionalType &functional)
      : space_( space ),
        functional_( functional )
    {
    }

    template <class DiscreteFunctionType>
    void assemble( DiscreteFunctionType &discreteFunction,
                   const int quadOrder = -1 ) const
    {
      const int quadratureOrder = ( quadOrder < 0 ) ? (2 * space_.order() + 1) : quadOrder;
      const DiscreteFunctionSpaceType &discreteFunctionSpace = space_;

      // set discreteFunction to zero
      discreteFunction.clear();

      Dune::TemporaryLocalFunction< DiscreteFunctionSpaceType >
        localFunctional( discreteFunctionSpace );

      // start grid traversal
      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin();
           it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        // init local function
        localFunctional.init( entity );

        // apply function locally
        functional_.applyLocal( localFunctional, quadratureOrder );

        // add to global DoF vector
        typedef typename DiscreteFunctionType :: LocalFunctionType  LocalFunctionType;
        LocalFunctionType dfLocalFunctional = discreteFunction.localFunction( entity );
        dfLocalFunctional += localFunctional;
      }

      // communicate data (for parallel runs)
      discreteFunction.communicate();
    }

    /*
    //! discreteFunction is an output parameter (kind of return value)
    template <class VectorType>
    void assemble( VectorType& functionalValues,
                   const int quadOrder = -1 ) const
    {
      // to be implemented for general vectors
      const int quadratureOrder = ( quadOrder < 0 ) ? (2 * space_.order() + 1) : quadOrder;
      const DiscreteFunctionSpaceType &discreteFunctionSpace = space_;

      // set discreteFunction to zero
      const size_t size = space_.size();
      for( size_t i = 0; i < size; ++i )
        functionalValues[ i ] = 0;

      Dune::TemporaryLocalFunction< DiscreteFunctionSpaceType >
        localFunctional( space_ );

      // start grid traversal
      const IteratorType endit = space_.end();
      for( IteratorType it = discreteFunctionSpace.begin();
           it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        // init local function
        localFunctional.init( entity );

        // apply function locally
        functional_.applyLocal( localFunctional, quadratureOrder );

        //const int localDofs = space.mapper().numEntityDofs(
        //for
        // add to global DoF vector
      }
    }
    */
  };

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_L2FUNCTIONAL_HH
