/**
  \file   main.cc
  \brief  Main file fir the finite element example.
  **/
// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

// dune-grid includes
#include <dune/grid/utility/gridtype.hh>

// dune-istl includes
#include <dune/istl/solvers.hh>

// dune-fem includes
//#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>

// reenable warnings
#include <dune/fem-tools/header/enablewarnings.hh>

// dune-functionals includes
#include <dune/functionals/container/factory.hh>
#include <dune/functionals/assembler/local/codim0/matrix.hh>
#include <dune/functionals/assembler/local/codim0/vector.hh>
#include <dune/functionals/assembler/system/affine.hh>
#include <dune/functionals/discretefunction/continuous.hh>
#include <dune/functionals/discretefunctionspace/continuous/lagrange.hh>
#include <dune/functionals/discretefunctionspace/subspace/linear.hh>
#include <dune/functionals/discretefunctionspace/subspace/affine.hh>
#include <dune/functionals/discreteoperator/local/codim0/integral.hh>
#include <dune/functionals/discretefunctional/local/codim0/integral.hh>

// dune-fem-functionals includes

// dune-fem-tools includes
#include <dune/fem-tools/common/string.hh>
#include <dune/fem-tools/common/printing.hh>
#include <dune/fem-tools/function/runtimefunction.hh>
//#include <dune/fem-tools/function/functiontools.hh>

using namespace Dune::Functionals;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

/**
  \brief  This represents the operation \f$fv\f$.

          \f$f\f$ is a given right hand side (in this case 1) and \f$v\f$ may be a local function, i.e. a
          testfunction.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$ and \f$v\f$ live in.
  **/
template< class FunctionSpaceImp >
class ProductEvaluation
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef Dune::FemTools::Function::Runtime< FunctionSpaceType >
    InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  ProductEvaluation( const std::string expression = "[1.0;1.0;1.0]" )
    : inducingFunction_( expression )
  {
  }

  //! copy constructor
  ProductEvaluation( const ProductEvaluation& other )
    : inducingFunction_( other.inducingFunction() )
  {
  }

  //! returns the inducing function
  const InducingFunctionType& inducingFunction() const
  {
    return inducingFunction_;
  }

  /**
    \brief      Evaluates \f$f(x)v(x)\f$ for a given local point \f$x\f$.

    \tparam     LocalTestFunctionType
                Type of the local function \f$v\f$, i.e. Dune::LocalFunction.
    \tparam     LocalPointType
                Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    \param[in]  localTestFunction
                The local function \f$v\f$.
    \param[in]  localPoint
                The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
                element.
    \return     \f$f(x)v(x)\f$
    **/
  template< class LocalTestBaseFunctionSetType, class LocalVectorType >
  void evaluate(  const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                  const DomainType& localPoint,
                  LocalVectorType& ret ) const
  {
    // get global point
    const DomainType globalPoint = localTestBaseFunctionSet.entity().geometry().global( localPoint );

    // evaluate inducing function
    RangeType functionValue( 0.0 );
    inducingFunction_.evaluate( globalPoint, functionValue );

    // evaluate set of local functions
    std::vector< RangeType > valuesLocalBaseFunctionSet( localTestBaseFunctionSet.size(), RangeType( 0.0 ) );
    localTestBaseFunctionSet.evaluate( localPoint, valuesLocalBaseFunctionSet );

    // do loop over all basis functions
    for( unsigned int i = 0; i < localTestBaseFunctionSet.size(); ++i )
    {
      ret[i] = functionValue * valuesLocalBaseFunctionSet[i];
    }
  }

private:

  const InducingFunctionType inducingFunction_;
}; // end class ProductEvaluation


/**
  \brief  This represents the operation \f$a\nabla u \nabla v\f$.

          \f$a\f$ is a given scalar function (in this case 1) and \f$u\f$ and \f$u\f$ may be local functions, i.e.
          ansatz- and testfunctions.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$, \f$u\f$ and \f$v\f$ live in.
  **/
template< class FunctionSpaceImp >
class EllipticEvaluation
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef Dune::FemTools::Function::Runtime< FunctionSpaceType >
    InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  EllipticEvaluation( const std::string expression = "[1.0;1.0;1.0]" )
    : inducingFunction_( expression )
  {
  }

  //! copy constructor
  EllipticEvaluation( const EllipticEvaluation& other )
    : inducingFunction_( other.inducingFunction() )
  {
  }

  //! returns the inducing function
  const InducingFunctionType& inducingFunction() const
  {
    return inducingFunction_;
  }

  /**
    * \brief      Evaluates \f$a(x)\nabla u(x) \nabla v(x)\f$ for a given local point \f$x\f$.
    *
    * \tparam     LocalAnsatzFunctionType
    *             Type of the local ansatz function \f$u\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalTestFunctionType
    *             Type of the local test function \f$v\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalPointType
    *             Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    * \param[in]  localAnsatzFunction
    *             The local function \f$u\f$.
    * \param[in]  localTestFunction
    *             The local function \f$v\f$.
    * \param[in]  localPoint
    *             The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
    *             element.
    * \return     \f$a(x)\nabla u(x) \nabla v(x)\f$
    **/
  template< class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType, class LocalMatrixType >
  void evaluate(  const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                  const LocalTestBaseFunctionSetType& localTestBaseFunctionSet,
                  const DomainType& localPoint,
                  LocalMatrixType& ret ) const
  {
//std::cout << "EllipticEvaluation::evaluate()" << std::endl;
//Dune::FemTools::Printing::printFieldVector( localPoint, "localPoint", std::cout );
    // get global point
    const DomainType globalPoint = localAnsatzBaseFunctionSet.entity().geometry().global( localPoint );
//Dune::FemTools::Printing::printFieldVector( globalPoint, "globalPoint", std::cout );

    // evaluate first gradient
    std::vector< JacobianRangeType > gradientLocalAnsatzBaseFunctionSet( localAnsatzBaseFunctionSet.size(), JacobianRangeType( 0.0 ) );
    localAnsatzBaseFunctionSet.jacobian( localPoint, gradientLocalAnsatzBaseFunctionSet );

    // evaluate second gradient
    std::vector< JacobianRangeType > gradientLocalTestBaseFunctionSet( localTestBaseFunctionSet.size(), JacobianRangeType( 0.0 ) );
    localTestBaseFunctionSet.jacobian( localPoint, gradientLocalTestBaseFunctionSet );

    // evaluate inducing function
    RangeType functionValue( 0.0 );
    inducingFunction_.evaluate( globalPoint, functionValue );

    for( unsigned int i = 0; i < localAnsatzBaseFunctionSet.size(); ++i )
    {
//std::cout << "i = " << i << " (of " << localAnsatzBaseFunctionSet.size() << ")" << std::endl;
//Dune::FemTools::Printing::printFieldVector( gradientLocalAnsatzBaseFunctionSet[i], "gradientLocalAnsatzBaseFunctionSet[" + Dune::FemTools::String::toString( i ) + "]", std::cout, "  " );
      for( unsigned int j = 0; j < localTestBaseFunctionSet.size(); ++j )
      {
//std::cout << "  j = " << j << " (of " << localTestBaseFunctionSet.size() << ")" << std::endl;
//Dune::FemTools::Printing::printFieldVector( gradientLocalTestBaseFunctionSet[i], "gradientLocalTestBaseFunctionSet[" + Dune::FemTools::String::toString( j ) + "]", std::cout, "    " );
        const RangeFieldType gradientProduct = gradientLocalAnsatzBaseFunctionSet[i][0] * gradientLocalTestBaseFunctionSet[j][0];
        ret[i][j] = functionValue * gradientProduct;
      }
    }
  }

private:

  const InducingFunctionType inducingFunction_;
}; // end class EllipticEvalaution


int main( int argc, char** argv )
{
  try{

    // MPI manager
    Dune::MPIManager::initialize( argc, argv );


    // grid
    static const unsigned int dimRange = 1;

    typedef Dune::GridSelector::GridType
      GridType;

    typedef Dune::AdaptiveLeafGridPart< GridType >
      GridPartType;

    const std::string dgfFilename = "../macrogrids/unitcube" + Dune::FemTools::String::toString( GRIDDIM ) + ".dgf";

    Dune::GridPtr< GridType > gridPtr( dgfFilename );

    GridPartType gridPart( *gridPtr );


    // function space
    typedef Dune::FunctionSpace< double, double, GridType::dimension, dimRange >
      FunctionSpaceType;

    typedef typename FunctionSpaceType::RangeFieldType
      RangeFieldType;


    // discrete function space
    typedef DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType, GridPartType, polOrder >
      DiscreteH1Type;

    const DiscreteH1Type discreteH1( gridPart );

    typedef DiscreteFunctionSpace::Subspace::Linear::Dirichlet< DiscreteH1Type >
      DiscreteH10Type;

    const DiscreteH10Type discreteH10( discreteH1 );

    typedef DiscreteFunctionSpace::Subspace::Affine::Dirichlet< DiscreteH10Type >
      DiscreteH1GType;

    const DiscreteH1GType discreteH1G( discreteH10, "[x+y;y;z]" );


    // local evaluation
    typedef ProductEvaluation< FunctionSpaceType >
      ProductEvaluationType;

    ProductEvaluationType productEvaluation( "[1.0;1.0;1.0]" );

    typedef EllipticEvaluation< FunctionSpaceType >
      EllipticEvaluationType;

    EllipticEvaluationType ellipticEvaluation( "[1.0;1.0;1.0]" );


    // operator and functional
    typedef DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
      LocalEllipticOperatorType;

    const LocalEllipticOperatorType localEllipticOperator( ellipticEvaluation );

    typedef DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
      LocalL2FunctionalType;

    const LocalL2FunctionalType localL2Functional( productEvaluation );

//    typedef typename LocalEllipticOperatorType::LocalFunctional< typename DiscreteH1GType::AffineShiftType >::Type
//      LocalAffineShiftFunctionalType;

//    const LocalAffineShiftFunctionalType localAffineShiftFunctional( localEllipticOperator, discreteH1G.affineShift() );

    // matrix, rhs and solution storage
    typedef Container::Matrix::Defaults< RangeFieldType, dimRange, dimRange >::BCRSMatrix
      MatrixFactory;

    typedef typename MatrixFactory::AutoPtrType
      MatrixPtrType;

    //! \todo the matrix factory should get two spaces (ansatz and test)
    MatrixPtrType A = MatrixFactory::create( discreteH1, discreteH1 );

    typedef Container::Vector::Defaults< RangeFieldType, dimRange >::BlockVector
      VectorFactory;

    typedef typename VectorFactory::AutoPtrType
      VectorPtrType;

    VectorPtrType F = VectorFactory::create( discreteH1 );

    VectorPtrType G = VectorFactory::create( discreteH1 );

    VectorPtrType u0 = VectorFactory::create( discreteH1 );


    // assembler
    typedef Assembler::Local::Codim0::Matrix< LocalEllipticOperatorType >
      LocalMatrixAssemblerType;

    const LocalMatrixAssemblerType localMatrixAssembler( localEllipticOperator );

    typedef Assembler::Local::Codim0::Vector< LocalL2FunctionalType >
      LocalVectorAssemblerType;

    const LocalVectorAssemblerType localVectorAssembler( localL2Functional );

    typedef Assembler::System::Affine< DiscreteH1GType, DiscreteH10Type >
      SystemAssemblerType;

    SystemAssemblerType systemAssembler( discreteH1G, discreteH10 );

    systemAssembler.assembleSystem( localMatrixAssembler, *A,
                                    localVectorAssembler, *F,
                                    *G );


    // preconditioner and solver
    typedef typename MatrixFactory::ContainerType
      MatrixContainerType;

    typedef typename VectorFactory::ContainerType
      VectorContainerType;

    typedef Dune::MatrixAdapter< MatrixContainerType, VectorContainerType, VectorContainerType >
      MatrixAdapterType;

    MatrixAdapterType matrix( *A );

    typedef Dune::SeqILU0< MatrixContainerType, VectorContainerType, VectorContainerType, 1 >
      PreconditionerType;

    PreconditionerType preconditioner( *A, 1.0 );

    typedef Dune::CGSolver< VectorContainerType >
      SolverType;

    SolverType solver( matrix, preconditioner, 1e-4, 100, 2 );

    Dune::InverseOperatorResult result;

    // u_0 = A^(-1) ( F - G )
    solver.apply( *u0, *F, result );


    // postprocessing
    typedef typename Dune::Functionals::DiscreteFunction::Continuous::BlockVector< DiscreteH1Type >
      DiscreteFunctionType;

    const DiscreteFunctionType solution( discreteH1, *u0, "solution" );
//    typedef Dune::Functionals::DiscreteFunctionSpace::Continuous::LagrangeFemAdapter< DiscreteH1Type >
//      LagrangeFemAdapterType;

//    const LagrangeFemAdapterType lagrangeFemAdapter( discreteH1 );

//    typedef Dune::BlockVectorDiscreteFunction< LagrangeFemAdapterType >
//      DiscreteFunctionType;

//    const DiscreteFunctionType solution( "solution", lagrangeFemAdapter );
//    DiscreteFunctionType solution = Dune::FemTools::Function::createFromVector( lagrangeFemAdapter, *u0 );

    Dune::FemTools::Function::writeToVTK( solution, "solution" );

    // done
    return 0;
  }
  catch( Dune::Exception& e )
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch( std::exception& e )
  {
    std::cerr << e.what() << std::endl;
  }
  catch( ... )
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
