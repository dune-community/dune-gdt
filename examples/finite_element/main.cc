#include <config.h>

#include <iostream>
#include <vector>

#include <dune/grid/utility/gridtype.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>

#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>


#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/constraints/dirichlet.hh>
#include <dune/fem/subspace/subspaces.hh>
#include <dune/fem/operator/ellipticfiniteelement.hh>
#include <dune/fem/functional/ltwo.hh>
#include <dune/fem/container/factory.hh>
#include <dune/fem/solver/femassembler.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

using namespace Dune::Functionals;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif


class GFunc
{
};

int main( int argc, const char *argv[] )
{

  typedef Dune::GridSelector::GridType
    HGridType;

  typedef Dune::AdaptiveLeafGridPart< HGridType >
    GridPartType;

  typedef Dune::FunctionSpace< double, double, HGridType::dimension, 1 >
    FunctionSpaceType;

  typedef Dune::Function< double, double >
    FunctionType;

  //typedef FunctionType
    ////InducingFunctionType;
  
  // function spaces
  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
    H1;

  typedef Dune::AdaptiveDiscreteFunction< H1 >
    DiscreteFunctionType;

  typedef DiscreteFunctionType
    InducingFunctionType;

  typedef Constraints::Dirichlet< H1 >
    DirichletConstraints;

  typedef Subspace::Linear< H1, DirichletConstraints >
    H10;

  typedef Operator::EllipticFiniteElement< H1, GFunc/*InducingFunctionType*/ >
    EllipticOperator;
  typedef Functional::L2< H1, InducingFunctionType >
    RHS;

  typedef Container::MatrixFactory< Dune::BCRSMatrix< double > >
    MatrixFactoryType;
  typedef MatrixFactoryType::AutoPtrType
    MatrixContainerPtr;
  typedef MatrixFactoryType::ContainerType
    MatrixContainer;

  typedef Container::VectorFactory< std::vector< double > >
    VectorFactoryType;
  typedef VectorFactoryType::AutoPtrType
    VectorContainerPtr;
  typedef VectorFactoryType::ContainerType
    VectorContainer;

  typedef Solver::FEMAssembler<MatrixContainer, VectorContainer>
    Assembler;
  typedef Dune::CGSolver< VectorContainer > 
    CG;
  typedef Dune::SeqILU0< MatrixContainer, VectorContainer, VectorContainer, 1 >
    SeqILU0Type;


  // FunctionType should be obsolete
  // It should be either a DoF-Container or an analytical function, and for
  // transition phase a discrete function.
  typedef Subspace::Affine< H10, GFunc >
    H1g;

  //create grid
  Dune::GridPtr< HGridType > gridPtr( "dummy.dgf" );

  //get grid part
  GridPartType gridPart( *gridPtr );

  //some functions
  GFunc                gFunc;
  GFunc                aFunc;
  GFunc                fFunc;

  //create spaces
  H1                   h1( gridPart );
  DirichletConstraints dirichletConstraints( h1 );
  H10                  h10( h1, dirichletConstraints );
  H1g                  h1g( h10, gFunc );


  EllipticOperator     ellipticFEM( h1, aFunc );
  RHS                  rhs( h1, fFunc );

  //SparsityPattern & pattern = h1.fullSparsityPattern();


  MatrixContainerPtr A  = MatrixFactoryType::create( h10 );
  VectorContainerPtr F  = VectorFactoryType::create( h10 );
  VectorContainerPtr G  = VectorFactoryType::create( h10 );

/*  MatrixContainer& A  = Container::MatrixFactory<MatrixContainer>::createRef( h1 );
 *  VectorContainer& F  = Container::VectorFactory<VectorContainer>::createRef( h1 );
 *  VectorContainer& G  = Container::VectorFactory<VectorContainer>::createRef( h1 );*/

  VectorContainerPtr u0 = VectorFactoryType::create( h10 );
  VectorContainerPtr u  = VectorFactoryType::create( h10 );


  Assembler::assembleMatrix( ellipticFEM, *A );
  Assembler::applyMatrixConstraints( h10, *A );
  Assembler::assembleVector( rhs, *F );
  Assembler::applyVectorConstraints( h10, *F );
  Assembler::assembleVector( ellipticFEM(gFunc), *G );


  CG cg(A, prec, 1e-10, 1000, true);

  InverseOperatorResult res;
  cg(u0, F-G, res);

  *u = *u0 + gFunc;

/*  DiscreteFunction dfU( h1, *u );

 *  dfU.evaluate( globalX );
 *  plot( dfU );*/

  return 0;
}
