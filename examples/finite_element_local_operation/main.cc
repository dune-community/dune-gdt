
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

// reenable warnings
#include <dune/fem-tools/header/enablewarnings.hh>

// dune-fem-functionals includes
#include <dune/fem/localoperation/interface.hh>
#include <dune/fem/localoperation/integrator.hh>
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


template< class FunctionSpaceImp >
class ProductOperation
  : public Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  template< class LocalFunctionType, class LocalPointType >
  RangeFieldType evaluate(  const LocalFunctionType& localFunction,
                            const LocalPointType& localPoint ) const
  {
    // init return value
    RangeFieldType ret = 0.0;

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
  * \brief  Represents the elliptic operation a(x) \gradient u(x) \gradient v(x) for given u, v, x.
  *         In this case, a = 1.
  **/
template< class FunctionSpaceImp >
class EllipticOperation
  : public Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
{
public:

  typedef FunctionSpaceImp
    FunctionSpaceType;

  typedef Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
    BaseType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  template< class AnsatzLocalFunctionType, class TestLocalFunctionType, class LocalPointType >
  RangeFieldType evaluate(  const AnsatzLocalFunctionType& ansatzLocalFunction,
                            const TestLocalFunctionType& testLocalFunction,
                            const LocalPointType& localPoint ) const
  {
    // init return value
    RangeFieldType ret = 0.0;

    // evaluate first gradient
    JacobianRangeType gradientAnsatzLocalFunction( 0.0 );
    ansatzLocalFunction.jacobian( localPoint, gradientAnsatzLocalFunction );

    // evaluate second gradient
    JacobianRangeType gradientTestLocalFunction( 0.0 );
    testLocalFunction.jacobian( localPoint, gradientTestLocalFunction );

    const RangeFieldType product = gradientAnsatzLocalFunction[0] * gradientTestLocalFunction[0];

    // 1.0 * \gradient u(x) \gradient v(x)
    ret = 1.0 * product;

    // return
    return ret;
  }

}; // end class EllipticOperation


// disable warnings about the dgf parser
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

    // local operations
    typedef ProductOperation< FunctionSpaceType >
      ProductOperationType;

    ProductOperationType productOperation;

    typedef EllipticOperation< FunctionSpaceType >
      EllipticOperationType;

    EllipticOperationType ellipticOperation;

    // integration
    typedef LocalOperation::VolumeIntegrator< FunctionSpaceType, ProductOperationType >
      ProductIntegratorType;

    ProductIntegratorType productIntegrator( productOperation );

    typedef LocalOperation::VolumeIntegrator< FunctionSpaceType, EllipticOperationType >
      EllipticIntegratorType;

    EllipticIntegratorType ellipticIntegrator( ellipticOperation );


    // discrete function space
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


    // operator and functional
    typedef Operator::FiniteElementLOP< DiscreteH1Type, DiscreteH1Type, EllipticIntegratorType >
      FEMellipticOperatorType;

    FEMellipticOperatorType femEllipticOperator( discreteH1, discreteH1, ellipticIntegrator );

    typedef Functional::FiniteElementLOP< DiscreteH1Type, ProductIntegratorType >
      FEMrhsFunctionalType;

    FEMrhsFunctionalType femRhsFunctional( discreteH1, productIntegrator );


    // matrix, rhs and solution storage
    typedef Dune::FieldMatrix< double, dimRange, dimRange >
      FieldMatrixType;

    typedef Container::MatrixFactory< Dune::BCRSMatrix< FieldMatrixType > >
      MatrixFactoryType;

    typedef MatrixFactoryType::ContainerType
      MatrixContainer;

    typedef MatrixFactoryType::AutoPtrType
      MatrixContainerPtr;

    typedef Container::VectorFactory< Dune::BlockVector< Dune::FieldVector< double, 1 > > >
      VectorFactoryType;

    typedef VectorFactoryType::ContainerType
      VectorContainer;

    typedef VectorFactoryType::AutoPtrType
      VectorContainerPtr;

    MatrixContainerPtr A  = MatrixFactoryType::create( discreteH10 );
    VectorContainerPtr F  = VectorFactoryType::create( discreteH10 );
    VectorContainerPtr u0 = VectorFactoryType::create( discreteH10 );


    // assembler
    typedef Assembler::FiniteElement< MatrixContainer, VectorContainer >
      Assembler;

    Assembler::assembleMatrix( femEllipticOperator, *A );
    Assembler::applyMatrixConstraints( discreteH10, *A );

    Assembler::assembleVector( femRhsFunctional, *F );
    Assembler::applyVectorConstraints( discreteH10, *F );


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
    cg.apply( *u0, *F, res );


    // postprocessing
    DiscreteFunctionType solution = Dune::FemTools::discreteFunctionFactory< DiscreteFunctionType >( discreteH1, *u0 );
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
