#include <config.h>

#include <iostream>
#include <vector>

#include <dune/grid/utility/gridtype.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/istl/bvector.hh>
#include <dune/common/fmatrix.hh>

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

/**
 * @class AFunc
 * function representing the coefficient a for the poisson problem
 * a \laplace u = f
 */
template <class FunctionSpaceImp>
class AFunc
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  void evaluate(const DomainType& x, RangeType& y) const
  {
    y = 1.0;
  }

  void jacobian(const DomainType& x, JacobianRangeType& y) const
  {
    y = 0.0;
  }
};

/**
 * @class FFunc
 * function representing f for the poisson problem
 * a \laplace u = f
 */
template <class FunctionSpaceImp>
class FFunc
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  void evaluate(const DomainType& x, RangeType& y) const
  {
    y = 0.0;
  }

  void jacobian(const DomainType& x, JacobianRangeType& y) const
  {
    y = 0.0;
  }
};

/**
 * @class GFunc
 * function representing the dirichlet data g for the poisson problem
 * a \laplace u = f
 */
template <class FunctionSpaceImp>
class GFunc
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  void evaluate(const DomainType& x, RangeType& y) const
  {
    y = 0.0;
  }

  void jacobian(const DomainType& x, JacobianRangeType& y) const
  {
    y = 0.0;
  }
};


int main(int argc, const char* argv[])
{
  static const unsigned int dimRange = 1;

  typedef Dune::GridSelector::GridType HGridType;

  typedef Dune::AdaptiveLeafGridPart<HGridType> GridPartType;

  typedef Dune::FunctionSpace<double, double, HGridType::dimension, dimRange> FunctionSpaceType;

  typedef Dune::Function<double, double> FunctionType;

  // typedef FunctionType
  ////InducingFunctionType;

  // function spaces
  typedef Dune::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, polOrder> H1;

  typedef Dune::AdaptiveDiscreteFunction<H1> DiscreteFunctionType;

  typedef DiscreteFunctionType InducingFunctionType;

  typedef Constraints::Dirichlet<H1> DirichletConstraints;

  typedef Subspace::Linear<H1, DirichletConstraints> H10;

  typedef Operator::EllipticFiniteElement<H1, AFunc<FunctionSpaceType> /*InducingFunctionType*/> EllipticOperator;
  typedef Functional::L2<H1, FFunc<FunctionSpaceType> /*InducingFunctionType*/> RHS;

  typedef Dune::FieldMatrix<double, dimRange, dimRange> FieldMatrixType;

  typedef Container::MatrixFactory<Dune::BCRSMatrix<FieldMatrixType>> MatrixFactoryType;
  typedef MatrixFactoryType::AutoPtrType MatrixContainerPtr;
  typedef MatrixFactoryType::ContainerType MatrixContainer;

  typedef Container::VectorFactory<Dune::BlockVector<Dune::FieldVector<double, 1>>> VectorFactoryType;
  typedef VectorFactoryType::AutoPtrType VectorContainerPtr;
  typedef VectorFactoryType::ContainerType VectorContainer;

  typedef Solver::FEMAssembler<MatrixContainer, VectorContainer> Assembler;
  typedef Dune::CGSolver<VectorContainer> CG;
  typedef Dune::SeqILU0<MatrixContainer, VectorContainer, VectorContainer, 1> SeqILU0Type;
  typedef Dune::MatrixAdapter<MatrixContainer, VectorContainer, VectorContainer> MatrixAdapterType;


  // FunctionType should be obsolete
  // It should be either a DoF-Container or an analytical function, and for
  // transition phase a discrete function.
  typedef Subspace::Affine<H10, GFunc<FunctionSpaceType>> H1g;

  // create grid
  Dune::GridPtr<HGridType> gridPtr("macrogrids/unitcube2.dgf");

  // get grid part
  GridPartType gridPart(*gridPtr);

  // some functions
  GFunc<FunctionSpaceType> gFunc;
  AFunc<FunctionSpaceType> aFunc;
  FFunc<FunctionSpaceType> fFunc;

  // create spaces
  H1 h1(gridPart);
  DirichletConstraints dirichletConstraints(h1);
  H10 h10(h1, dirichletConstraints);
  H1g h1g(h10, gFunc);


  EllipticOperator ellipticFEM(h1, aFunc);
  RHS rhs(h1, fFunc);

  // SparsityPattern & pattern = h1.fullSparsityPattern();


  MatrixContainerPtr A = MatrixFactoryType::create(h10);
  VectorContainerPtr F = VectorFactoryType::create(h10);
  VectorContainerPtr G = VectorFactoryType::create(h10);

  /*  MatrixContainer& A  = Container::MatrixFactory<MatrixContainer>::createRef( h1 );
   *  VectorContainer& F  = Container::VectorFactory<VectorContainer>::createRef( h1 );
   *  VectorContainer& G  = Container::VectorFactory<VectorContainer>::createRef( h1 );*/

  VectorContainerPtr u0 = VectorFactoryType::create(h10);
  VectorContainerPtr u  = VectorFactoryType::create(h10);


  Assembler::assembleMatrix(ellipticFEM, *A);
  Assembler::applyMatrixConstraints(h10, *A);
  Assembler::assembleVector(rhs, *F);
  Assembler::applyVectorConstraints(h10, *F);
  Assembler::assembleVector(ellipticFEM(gFunc), *G);

  MatrixAdapterType op(*A);
  SeqILU0Type prec(*A, 1.0);

  CG cg(op, prec, 1e-10, 1000, 1);

  Dune::InverseOperatorResult res;
  *F -= *G;
  cg.apply(*u0, *F, res);

  // @todo implement gFunc
  //*u = *u0 + gFunc;

  /*  DiscreteFunction dfU( h1, *u );

   *  dfU.evaluate( globalX );
   *  plot( dfU );*/

  return 0;
}
