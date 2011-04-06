
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>

// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

// dune-grid includes
#include <dune/grid/utility/gridtype.hh>

// dune-istl includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/bvector.hh>

// dune-fem includes
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// reenable warnings about problems in dune headers
#include <dune/fem-tools/header/enablewarnings.hh>

// dune-fem-functionals includes
#include <dune/fem/constraints/dirichlet.hh>
#include <dune/fem/subspace/subspaces.hh>
#include <dune/fem/operator/finiteelement.hh>
#include <dune/fem/functional/finiteelement.hh>
#include <dune/fem/container/factory.hh>
#include <dune/fem/solver/femassembler.hh>

// dune-fem-tools includes
#include <dune/fem-tools/function/functiontools.hh>
#include <dune/fem-tools/space/projection.hh>

using namespace Dune::Functionals;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

/**
  * \brief  Represents the elliptic operation a(x) \gradient u(x) \gradient v(x) for given u, v, x.
  *         In this case, a = 1.
  **/
class EllipticOperation
{
public:

  template< class FirstLocalFunctionType, class SecondLocalFunctionType, class LocalPointType >
  double operate( const FirstLocalFunctionType& firstLocalFunction,
                  const SecondLocalFunctionType& secondLocalFunction,
                  const LocalPointType& localPoint ) const
  {
    // init return value
    double ret = 0.0;

    // some types we will need
    typedef typename SecondLocalFunctionType::EntityType
      EntityType;

    typedef typename FirstLocalFunctionType::JacobianRangeType
      JacobianRangeType;

    // first gradient
    JacobianRangeType firstGradient( 0.0 );
    firstLocalFunction.jacobian( localPoint, firstGradient );

    // second gradient
    JacobianRangeType secondGradient( 0.0 );
    secondLocalFunction.jacobian( localPoint, secondGradient );

    const double product = firstGradient[0] * secondGradient[0];

    // 1.0 * \gradient u(x) \gradient v(x)
    ret = 1.0 * product;

    // return
    return ret;
  }

}; // end class EllipticOperation

/**
  * \brief  Represents the product operation f(x) v(x) for given v, x.
  *         In this case, f = 1.
  **/
class ProductOperation
{
public:

  template< class LocalFunctionType, class LocalPointType >
  double operate( const LocalFunctionType& localFunction,
                  const LocalPointType& localPoint ) const
  {
    // init return value
    double ret = 0.0;

    // some types we will need
    typedef typename LocalFunctionType::RangeType
      RangeType;

    // evaluate local function
    RangeType localFunctionEvaluated( 0.0 );
    localFunction.evaluate( localPoint, localFunctionEvaluated );

    // 1.0 * v(x)
    ret = 1.0 * localFunctionEvaluated;

    // return
    return ret;
  }

}; // end class ProductOperation

/**
 * @class AnalyticalDirichletValues
 * function representing the dirichlet data g for the poisson problem
 * a \laplace u = f
 */
template< class FunctionSpaceImp >
class AnalyticalDirichletValues
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  void evaluate( const DomainType& arg, RangeType& ret ) const
  {
    ret = 1.0;

    if( !(arg[0] > 0.0) )
    {
      // bottom
      ret = 1.0;
    }
    else if( !(arg[0] < 1.0) )
    {
      // top
      ret = 1.0;
    }
    if( !(arg[1] > 0.0) )
    {
      // left
      ret = 1.0;
    }
    else if( !(arg[1] < 1.0) )
    {
      // right
      ret = 1.0;
    }
  }

  void jacobian( const DomainType& arg, JacobianRangeType& ret ) const
  {
    ret = 0.0;
  }

};


// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>


int main( int argc, char** argv )
{
  try{

    // MPI manager
    Dune::MPIManager::initialize ( argc, argv );


    // grid
    static const unsigned int dimRange = 1;

    typedef Dune::GridSelector::GridType
      HGridType;

    typedef Dune::AdaptiveLeafGridPart< HGridType >
      GridPartType;

    Dune::GridPtr< HGridType > gridPtr( "macrogrids/unitcube2.dgf" );

    GridPartType gridPart( *gridPtr );


    // function spaces and functions
    typedef Dune::FunctionSpace< double, double, HGridType::dimension, dimRange >
      FunctionSpaceType;

    typedef Dune::Function< double, double >
      FunctionType;

    typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteH1Type;

    DiscreteH1Type discreteH1( gridPart );

    typedef Dune::AdaptiveDiscreteFunction< DiscreteH1Type >
      DiscreteFunctionType;

    typedef Constraints::Dirichlet< DiscreteH1Type >
      DirichletConstraints;

    DirichletConstraints dirichletConstraints( discreteH1 );

    typedef Subspace::Linear< DiscreteH1Type, DirichletConstraints >
      DiscreteH10Type;

    DiscreteH10Type discreteH10( discreteH1, dirichletConstraints );

    typedef AnalyticalDirichletValues< FunctionSpaceType >
      AnalyticalDirichletValuesType;

    AnalyticalDirichletValuesType analyticalDirichletValues;

    DiscreteFunctionType discreteDirichletValues( "dirichlet_values", discreteH1 );

    Dune::FemTools::BetterL2Projection::project( analyticalDirichletValues, discreteDirichletValues );

    // FunctionType should be obsolete
    // It should be either a DoF-Container or an analytical function, and for
    // transition phase a discrete function.
    typedef Subspace::Affine< DiscreteH10Type, DiscreteFunctionType >
      DiscreteH1gType;

    DiscreteH1gType discreteH1g( discreteH10, discreteDirichletValues );


    // operator and functional
    typedef Operator::FiniteElement< DiscreteH1Type, EllipticOperation >
      FEMellipticOperator;

    EllipticOperation ellipticOperation;

    FEMellipticOperator femEllipticOperator( discreteH1, ellipticOperation );

    typedef Functional::FiniteElement< DiscreteH1Type, ProductOperation >
      FEMrhsFunctional;

    ProductOperation productOperation;

    FEMrhsFunctional femRhsFunctional( discreteH1, productOperation );


    // matrix, rhs and solution storage
    typedef Dune::FieldMatrix< double, dimRange, dimRange >
      FieldMatrixType;

    typedef Container::MatrixFactory< Dune::BCRSMatrix< FieldMatrixType > >
      MatrixFactoryType;

    typedef MatrixFactoryType::AutoPtrType
      MatrixContainerPtr;

    typedef MatrixFactoryType::ContainerType
      MatrixContainer;

    typedef Container::VectorFactory< Dune::BlockVector< Dune::FieldVector< double, 1 > > >
      VectorFactoryType;

    typedef VectorFactoryType::AutoPtrType
      VectorContainerPtr;

    typedef VectorFactoryType::ContainerType
      VectorContainer;

//    SparsityPattern & pattern = h1.fullSparsityPattern();

    MatrixContainerPtr A  = MatrixFactoryType::create( discreteH10 );
    VectorContainerPtr F  = VectorFactoryType::create( discreteH10 );
    VectorContainerPtr G  = VectorFactoryType::create( discreteH1 );

//    MatrixContainer& A  = Container::MatrixFactory<MatrixContainer>::createRef( discreteH10 );
//    VectorContainer& F  = Container::VectorFactory<VectorContainer>::createRef( discreteH10 );
//    VectorContainer& G  = Container::VectorFactory<VectorContainer>::createRef( discreteH1 );

    VectorContainerPtr u0 = VectorFactoryType::create( discreteH10 );


    // assembler
    typedef Solver::FEMAssembler< MatrixContainer, VectorContainer >
      Assembler;

    Assembler::assembleMatrix( femEllipticOperator, *A );
    Assembler::applyMatrixConstraints( discreteH10, *A );

    Assembler::assembleVector( femRhsFunctional, *F );
    Assembler::applyVectorConstraints( discreteH10, *F );

//    Assembler::assembleVector( femEllipticOperator( discreteDirichletValues ), *G );


    // preconditioner and solver
    typedef Dune::MatrixAdapter< MatrixContainer, VectorContainer, VectorContainer >
      MatrixAdapterType;

    MatrixAdapterType op( *A );

    typedef Dune::SeqILU0< MatrixContainer, VectorContainer, VectorContainer, 1 >
      SeqILU0Type;

    SeqILU0Type prec( *A, 1.0 );

    typedef Dune::CGSolver< VectorContainer >
      CG;

    CG cg( op, prec, 1e-4, 100, 2 );

    Dune::InverseOperatorResult res;

    // u_0 = A^(-1) ( F - G )
//    *G = G->operator*=( -1.0 );
//    *F = F->operator+=( *G );
    cg.apply( *u0, *F, res );


    // postprocessing
    DiscreteFunctionType solution = Dune::FemTools::discreteFunctionFactory< DiscreteFunctionType >( discreteH1, *u0 );
//    solution += discreteDirichletValues;
    Dune::FemTools::writeDiscreteFunctionToVTK( solution, "solution" );

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
