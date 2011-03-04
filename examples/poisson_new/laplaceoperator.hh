#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

//- Dune includes
#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem-howto/probleminterfaces.hh>

namespace Dune
{

  //! \brief The Laplace operator
  template< class DiscreteFunction, class MatrixTraits >             /*@LST0S@*/
  class LaplaceOperator
  : public Operator< typename DiscreteFunction::RangeFieldType,
                     typename DiscreteFunction::RangeFieldType,
                     DiscreteFunction,
                     DiscreteFunction >,
    public OEMSolver::PreconditionInterface                         /*@\label{poi:precif}@*/
  {                                                                 /*@LST0E@*/
    typedef LaplaceOperator< DiscreteFunction, MatrixTraits > ThisType;
    typedef Operator< typename DiscreteFunction::RangeFieldType, typename DiscreteFunction::RangeFieldType,
                      DiscreteFunction, DiscreteFunction > BaseType;

    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    
  public:
    //! type of discrete functions
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of problem 
    typedef ProblemInterface< typename DiscreteFunctionSpaceType :: FunctionSpaceType > 
        ProblemType; 

    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;
       
  protected:
    //! type of jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;
    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  public:
    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };
        
    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

    typedef typename MatrixTraits
      :: template  MatrixObject< LagrangeMatrixTraits< MatrixTraits > >
      :: MatrixObjectType LinearOperatorType;

    //! get important types from the MatrixObject 
    typedef typename LinearOperatorType :: LocalMatrixType LocalMatrixType;
    typedef typename LinearOperatorType :: PreconditionMatrixType PreconditionMatrixType;
    typedef typename LinearOperatorType :: MatrixType MatrixType;

  protected:
    // type of DofManager
    typedef DofManager< GridType > DofManagerType;
   
  public:
    //! constructor 
    LaplaceOperator ( const DiscreteFunctionSpaceType &dfSpace, const ProblemType &problem )
    : discreteFunctionSpace_( dfSpace ),
      dofManager_( DofManagerType::instance( dfSpace.grid() ) ),
      linearOperator_( discreteFunctionSpace_, discreteFunctionSpace_ ),
      sequence_( -1 ),
      problem_( problem ),
      gradCache_( discreteFunctionSpace_.mapper().maxNumDofs() ),
      gradDiffusion_( discreteFunctionSpace_.mapper().maxNumDofs() )
    {}
        
  private:
    // prohibit copying
    LaplaceOperator ( const ThisType & );

  public:                                                           /*@LST0S@*/
    //! apply the operator
    virtual void operator() ( const DiscreteFunctionType &u, 
                              DiscreteFunctionType &w ) const 
    {
      systemMatrix().apply( u, w );                                 /*@\label{poi:matrixEval}@*/
    }                                                               /*@LST0E@*/
  
    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType &preconditionMatrix () const
    {
      return systemMatrix().preconditionMatrix();
    }

    //! return true if preconditioning is enabled
    bool hasPreconditionMatrix () const
    {
      return linearOperator_.hasPreconditionMatrix();
    }

    //! print the system matrix into a stream
    void print ( std :: ostream & out = std :: cout ) const 
    {
      systemMatrix().matrix().print( out );
    }

    //! return reference to discreteFunctionSpace
    const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    //! return reference to problem
    const ProblemType& problem() const
    {
      return problem_;
    } 

    /*! \brief obtain a reference to the system matrix
     *
     *  The assembled matrix is returned. If the system matrix has not been
     *  assembled, yet, the assembly is performed.
     *
     *  \returns a reference to the system matrix
     */
    LinearOperatorType &systemMatrix () const                        /*@LST0S@*/
    {
      // if stored sequence number it not equal to the one of the
      // dofManager (or space) then the grid has been changed
      // and matrix has to be assembled new
      if( sequence_ != dofManager_.sequence() )                     /*@\label{poi:sequence}@*/
        assemble();

      return linearOperator_;
    }

    /** \brief perform a grid walkthrough and assemble the global matrix */
    void assemble () const 
    {
      const DiscreteFunctionSpaceType &space = discreteFunctionSpace();

      // reserve memory for matrix 
      linearOperator_.reserve();                                   /*@LST0E@*/

      // create timer (also stops time)
      Timer timer;

      // clear matrix                                             /*@LST0S@*/ 
      linearOperator_.clear();

      // apply local matrix assembler on each element
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      IteratorType end = space.end();
      for(IteratorType it = space.begin(); it != end; ++it)
      {
        assembleLocal( *it );                                       /*@\label{poi:applAss}@*/
      }

      // get elapsed time 
      const double assemblyTime = timer.elapsed();
      // in verbose mode print times 
      if ( Parameter :: verbose () )
        std :: cout << "Time to assemble matrix: " << assemblyTime << "s" << std :: endl;

      // get grid sequence number from space (for adaptive runs)    /*@LST0S@*/
      sequence_ = dofManager_.sequence();
    }

  protected:
    //! assemble local matrix for given entity
    template< class EntityType >
    void assembleLocal( const EntityType &entity ) const
    {
      // extract type of geometry from entity 
      typedef typename EntityType :: Geometry Geometry;

      // assert that matrix is not build on ghost elements 
      assert( entity.partitionType() != GhostEntity );

      // cache geometry of entity 
      const Geometry &geometry = entity.geometry();

      // get local matrix from matrix object
      LocalMatrixType localMatrix                                   /*@\label{poi:localMInit}@*/
        = linearOperator_.localMatrix( entity, entity );
      
      // get base function set 
      const BaseFunctionSetType &baseSet = localMatrix.domainBaseFunctionSet();/*@\label{poi:baseSetInit}@*/

      // get number of local base functions 
      const size_t numBaseFunctions = baseSet.numBaseFunctions();
            
      // create quadrature of appropriate order 
      QuadratureType quadrature( entity, 2 * (polynomialOrder - 1) );/*@\label{poi:quadraInit}@*/

      // loop over all quadrature points
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // get local coordinate of quadrature point 
        const typename QuadratureType :: CoordinateType &x
          = quadrature.point( pt );

        // get jacobian inverse transposed 
        const typename Geometry :: Jacobian& inv
          = geometry.jacobianInverseTransposed( x );

        // extract type of diffusion coefficient from problem
        typedef typename ProblemType :: DiffusionMatrixType DiffusionMatrixType;
        DiffusionMatrixType K;

        // evaluate diffusion matrix 
        problem().K( geometry.global( x ), K );

        // for all base functions evaluate the gradient
        // on quadrature point pt and apply jacobian inverse 
        baseSet.jacobianAll( quadrature[ pt ], inv, gradCache_ );

        // apply diffusion tensor  
        for( size_t i = 0; i < numBaseFunctions; ++i )
          K.mv( gradCache_[ i ][ 0 ], gradDiffusion_[ i ][ 0 ] );
        
        // evaluate integration weight 
        weight_ = quadrature.weight( pt ) * geometry.integrationElement( x );
        
        // add scalar product of gradients to local matrix 
        updateLocalMatrix( localMatrix );
      }
    }
  
    //! add scalar product of cached gradients to local matrix
    void updateLocalMatrix ( LocalMatrixType &localMatrix ) const
    {
      const size_t rows    = localMatrix.rows();
      const size_t columns = localMatrix.columns();
      for( size_t i = 0; i < rows; ++i )
      {
        for ( size_t j = 0; j < columns; ++j )
        {
          const RangeFieldType value
            = weight_ * (gradCache_[ i ][ 0 ] * gradDiffusion_[ j ][ 0 ]);
          localMatrix.add( j, i, value );
        }
      }
    }                                                            /*@LST0E@*/
  
  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const DofManagerType &dofManager_;

    //! pointer to the system matrix
    mutable LinearOperatorType linearOperator_;
 
    //! flag indicating whether the system matrix has been assembled
    mutable int sequence_;
      
    //! the diffusion tensor 
    const ProblemType& problem_;

    mutable DynamicArray< JacobianRangeType > gradCache_;
    mutable DynamicArray< JacobianRangeType > gradDiffusion_;
    mutable RangeFieldType weight_;    
  };                                                                  /*@LST0S@*//*@LST0E@*/

} // end namespace 
#endif
