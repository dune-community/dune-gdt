#include <iostream>

#include <dune/grid/utility/gridtype.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/constraints/dirichlet.hh>
#include <dune/fem/subspace/subspaces.hh>
#include <dune/fem/operator/ellipticfiniteelement.hh>
#include <dune/fem/functional/ltwo.hh>
#include <dune/fem/container/factory.hh>
#include <dune/fem/solver/femassembler.hh>

using namespace Dune::Functionals;

class gFunc
{
};

int main(int argc, const char *argv[])
{

  typedef Dune::AdaptiveLeafGridPart<Dune::GridSelector::GridType>
    GridPartType;

  // function spaces
  typedef Dune::LagrangeDiscreteFunctionSpace<FunctionSpace, GridPartType, polOrder>
    H1;

  typedef Constraints::Dirichlet<H1>
    DirichletConstraints;

  typedef Subspace::Linear<H1, DirichletConstraints>
    H10;

  // FunctionType should be obsolete
  // It should be either a DoF-Container or an analytical function, and for
  // transition phase a discrete function.
  typedef Subspace::Affine<H10, FunctionType>
    H1g;

  H1                   h1( gridPart );
  DirichletConstraints dirichletConstraints( h1 );
  H10                  h10( h1, dirichletConstraints );
  H1g                  h1g( h10, gFunc );


  typedef Operator::EllipticFiniteElement< H1, InducingFunction >
    EllipticFEM;

  typedef Functional::L2< H1, InducingFunction >
    RHS;

  EllipticFEM ellipticFEM( h1, aFunc );
  RHS         rhs( h1, fFunc );

  //SparsityPattern & pattern = h1.fullSparsityPattern();

  typedef Container::MatrixFactory<Dune::BCRSMatrix<double> >
    MatrixFactoryType;
  typedef typename MatrixFactoryType::AutoPtrType
    MatrixContainerPtr;
  typedef typename MatrixFactoryType::ContainerType
    MatrixContainer;

  typedef Container::VectorFactory<Dune::BCRSVector<double> >
    VectorFactoryType;
  typedef typename VectorFactoryType::AutoPtrType
    VectorContainerPtr;
  typedef typename VectorFactoryType::ContainerType
    VectorContainer;

  MatrixContainerPtr A  = MatrixFactoryType::create( h10 );
  VectorContainerPtr F  = VectorFactoryType::create( h10 );
  VectorContainerPtr G  = VectorFactoryType::create( h10 );

/*  MatrixContainer& A  = Container::MatrixFactory<MatrixContainer>::createRef( h1 );
 *  VectorContainer& F  = Container::VectorFactory<VectorContainer>::createRef( h1 );
 *  VectorContainer& G  = Container::VectorFactory<VectorContainer>::createRef( h1 );*/

  VectorContainerPtr u0 = VectorFactoryType::create( h10 );
  VectorContainerPtr u  = VectorFactoryType::create( h10 );

  typedef Solver::FEMAssembler<MatrixContainer, VectorContainer>
    Assembler;

  Assembler::assembleMatrix( ellipticFEM, *A );
  Assembler::applyMatrixConstraints( h10, *A );
  Assembler::assembleVector( rhs, *F );
  Assembler::applyVectorConstraints( h10, *F );
  Assembler::assembleVector( ellipticFEM(gFunc), *G );
  Assembler::applyVectorConstraints( h10, *G );

  typedef Dune::CGSolver< VectorContainer > CG;
  typedef Dune::SeqILU0< MatrixContainer, VectorContainer, VectorContainer, int >

  CG cg(A, prec, 1e-10, 1000, true);

  InverseOperatorResult res;
  cg(u0, F-G, res);

  *u = *u0 + gFunc;

/*  DiscreteFunction dfU( h1, *u );

 *  dfU.evaluate( globalX );
 *  plot( dfU );*/

  return 0;
}
