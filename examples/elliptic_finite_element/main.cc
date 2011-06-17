/**
  \file   main.cc
  \brief  Main file fir the finite element example.
  **/
// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>

#ifdef HAVE_CONFIG_H
#include "config.h"
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
#include <dune/fem/function/adaptivefunction.hh>
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
#include <dune/functionals/discretefunctionspace/continuous/lagrange.hh>
#include <dune/functionals/discretefunctionspace/subspace/linear.hh>
#include <dune/functionals/discretefunctionspace/subspace/affine.hh>
#include <dune/functionals/discreteoperator/local/codim0/integral.hh>
#include <dune/functionals/discretefunctional/local/codim0/integral.hh>

// dune-fem-tools includes
#include <dune/fem-tools/common/string.hh>
#include <dune/fem-tools/function/runtimefunction.hh>
#include <dune/fem-tools/function/functiontools.hh>

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
template <class FunctionSpaceImp>
class ProductEvaluation
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  ProductEvaluation(const std::string expression = "[1.0;1.0;1.0]")
    : inducingFunction_(expression)
  {
  }

  //! copy constructor
  ProductEvaluation(const ProductEvaluation& other)
    : inducingFunction_(other.inducingFunction())
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
  template <class LocalTestFunctionType>
  const RangeFieldType evaluate(const LocalTestFunctionType& localTestFunction, const DomainType& localPoint) const
  {
    // get global point
    const DomainType globalPoint = localTestFunction.entity().geometry().global(localPoint);

    // evaluate local function
    RangeType localTestFunctionValue(0.0);
    localTestFunction.evaluate(localPoint, localTestFunctionValue);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_.evaluate(globalPoint, functionValue);

    // f(x) * v(x)
    const RangeFieldType ret = functionValue * localTestFunctionValue;

    return ret;
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
template <class FunctionSpaceImp>
class EllipticEvaluation
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> InducingFunctionType;

  //! constructor, takes the inducing functions expression as a runtime parameter
  EllipticEvaluation(const std::string expression = "[1.0;1.0;1.0]")
    : inducingFunction_(expression)
  {
  }

  //! copy constructor
  EllipticEvaluation(const EllipticEvaluation& other)
    : inducingFunction_(other.inducingFunction())
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
  template <class LocalAnsatzFunctionType, class LocalTestFunctionType>
  const RangeFieldType evaluate(const LocalAnsatzFunctionType& localAnsatzFunction,
                                const LocalTestFunctionType& localTestFunction, const DomainType& localPoint) const
  {
    // get global point
    const DomainType globalPoint = localAnsatzFunction.entity().geometry().global(localPoint);

    // evaluate first gradient
    JacobianRangeType gradientLocalAnsatzFunction(0.0);
    localAnsatzFunction.jacobian(localPoint, gradientLocalAnsatzFunction);

    // evaluate second gradient
    JacobianRangeType gradientLocalTestFunction(0.0);
    localTestFunction.jacobian(localPoint, gradientLocalTestFunction);

    const RangeType gradientProduct = gradientLocalAnsatzFunction[0] * gradientLocalTestFunction[0];

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_.evaluate(globalPoint, functionValue);

    // a(x) * \gradient u(x) * \gradient v(x)
    const RangeFieldType ret = functionValue * gradientProduct;

    return ret;
  }

private:
  const InducingFunctionType inducingFunction_;
}; // end class EllipticEvalaution


int main(int argc, char** argv)
{
  try {

    // MPI manager
    Dune::MPIManager::initialize(argc, argv);


    // grid
    static const unsigned int dimRange = 1;

    typedef Dune::GridSelector::GridType GridType;

    typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;

    const std::string dgfFilename = "../macrogrids/unitcube" + Dune::FemTools::String::toString(GRIDDIM) + ".dgf";

    Dune::GridPtr<GridType> gridPtr(dgfFilename);

    GridPartType gridPart(*gridPtr);


    // function space
    typedef Dune::FunctionSpace<double, double, GridType::dimension, dimRange> FunctionSpaceType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;


    // discrete function space
    typedef DiscreteFunctionSpace::Continuous::Lagrange<FunctionSpaceType, GridPartType, polOrder> DiscreteH1Type;

    const DiscreteH1Type discreteH1(gridPart);

    typedef DiscreteFunctionSpace::Subspace::Linear::Dirichlet<DiscreteH1Type> DiscreteH10Type;

    const DiscreteH10Type discreteH10(discreteH1);

    typedef DiscreteFunctionSpace::Subspace::Affine::Dirichlet<DiscreteH10Type> DiscreteH1GType;

    const DiscreteH1GType discreteH1G(discreteH10, "[x+y;y;z]");


    // local evaluation
    typedef ProductEvaluation<FunctionSpaceType> ProductEvaluationType;

    ProductEvaluationType productEvaluation("[1.0;1.0;1.0]");

    typedef EllipticEvaluation<FunctionSpaceType> EllipticEvaluationType;

    EllipticEvaluationType ellipticEvaluation("[1.0;1.0;1.0]");


    // operator and functional
    typedef DiscreteOperator::Local::Codim0::Integral<EllipticEvaluationType> LocalEllipticOperatorType;

    const LocalEllipticOperatorType localEllipticOperator(ellipticEvaluation);

    typedef DiscreteFunctional::Local::Codim0::Integral<ProductEvaluationType> LocalL2FunctionalType;

    const LocalL2FunctionalType localL2Functional(productEvaluation);

    typedef typename LocalEllipticOperatorType::LocalFunctional<typename DiscreteH1GType::AffineShiftType>::Type
        LocalAffineShiftFunctionalType;

    const LocalAffineShiftFunctionalType localAffineShiftFunctional(localEllipticOperator, discreteH1G.affineShift());

    // matrix, rhs and solution storage
    //! \todo the matrix factory should get two spaces (ansatz and test)
    typedef Container::Matrix::Defaults<RangeFieldType, dimRange, dimRange>::BCRSMatrix MatrixFactory;

    typedef typename MatrixFactory::AutoPtrType MatrixPtrType;

    MatrixPtrType A = MatrixFactory::create(discreteH1);

    typedef Container::Vector::Defaults<RangeFieldType, dimRange>::BlockVector VectorFactory;

    typedef typename VectorFactory::AutoPtrType VectorPtrType;

    VectorPtrType F = VectorFactory::create(discreteH1);

    VectorPtrType G = VectorFactory::create(discreteH1);

    VectorPtrType u0 = VectorFactory::create(discreteH1);


    // assembler
    typedef Assembler::Local::Codim0::Matrix<LocalEllipticOperatorType> LocalMatrixAssemblerType;

    const LocalMatrixAssemblerType localMatrixAssembler(localEllipticOperator);

    typedef Assembler::Local::Codim0::Vector<LocalL2FunctionalType> LocalVectorAssemblerType;

    const LocalVectorAssemblerType localVectorAssembler(localL2Functional);

    typedef Assembler::System::Affine<DiscreteH1GType, DiscreteH10Type> SystemAssemblerType;

    SystemAssemblerType systemAssembler(discreteH1G, discreteH10);

    systemAssembler.assembleSystem(localMatrixAssembler, *A, localVectorAssembler, *F, *G);


    // preconditioner and solver
    typedef typename MatrixFactory::ContainerType MatrixContainerType;

    typedef typename VectorFactory::ContainerType VectorContainerType;

    typedef Dune::MatrixAdapter<MatrixContainerType, VectorContainerType, VectorContainerType> MatrixAdapterType;

    MatrixAdapterType matrix(*A);

    typedef Dune::SeqILU0<MatrixContainerType, VectorContainerType, VectorContainerType, 1> PreconditionerType;

    PreconditionerType preconditioner(*A, 1.0);

    typedef Dune::CGSolver<VectorContainerType> SolverType;

    SolverType solver(matrix, preconditioner, 1e-4, 100, 2);

    Dune::InverseOperatorResult result;

    // u_0 = A^(-1) ( F - G )
    solver.apply(*u0, *F, result);


    // postprocessing
    typedef Dune::AdaptiveDiscreteFunction<typename DiscreteH1Type::HostSpaceType> DiscreteFunctionType;

    DiscreteFunctionType solution =
        Dune::FemTools::Function::createFromVector<DiscreteFunctionType>(discreteH1.hostSpace(), *u0);

    Dune::FemTools::Function::writeToVTK(solution, "solution");

    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
