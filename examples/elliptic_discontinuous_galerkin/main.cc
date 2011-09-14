/**
  \file   main.cc
  \brief  Main file for the elliptic discontinuous galerkin example.
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
//#include <dune/common/fmatrix.hh>

// dune-grid includes
#include <dune/grid/utility/gridtype.hh>

// dune-istl includes
#include <dune/istl/solvers.hh>

// dune-fem includes
//#include <dune/fem/function/adaptivefunction.hh>
//#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/gridpart/gridpart.hh>
//#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>

// reenable warnings
#include <dune/fem-tools/header/enablewarnings.hh>

// dune-functionals includes
#include <dune/functionals/discretefunctionspace/discontinuous/orthonormal.hh>
//#include <dune/functionals/container/factory.hh>
//#include <dune/functionals/assembler/local/codim0/matrix.hh>
//#include <dune/functionals/assembler/local/codim0/vector.hh>
//#include <dune/functionals/assembler/system/unconstrained.hh>
//#include <dune/functionals/discretefunction/continuous.hh>
//#include <dune/functionals/discretefunction/femadapter.hh>
//#include <dune/functionals/discreteoperator/local/codim0/integral.hh>
//#include <dune/functionals/discretefunctional/local/codim0/integral.hh>
//#include <dune/functionals/evaluation/unary/scale.hh>
//#include <dune/functionals/evaluation/binary/elliptic.hh>

// dune-fem-tools includes
#include <dune/fem-tools/common/string.hh>
#include <dune/fem-tools/common/printing.hh>
//#include <dune/fem-tools/function/runtimefunction.hh>
#include <dune/fem-tools/discretefunction.hh>

using namespace Dune::Functionals;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif


int main(int argc, char** argv)
{
  try {

    // MPI manager
    Dune::MPIManager::initialize(argc, argv);


    // grid
    typedef Dune::GridSelector::GridType GridType;

    typedef Dune::LeafGridPart<GridType> GridPartType;

    const std::string dgfFilename = "../macrogrids/unitcube" + Dune::FemTools::String::toString(GRIDDIM) + ".dgf";

    Dune::GridPtr<GridType> gridPtr(dgfFilename);

    GridPartType gridPart(*gridPtr);

    static const unsigned int dimDomain = GridType::dimension;

    static const unsigned int dimRange = 1;

    // function space
    typedef Dune::FunctionSpace<double, double, dimDomain, dimRange> FunctionSpaceType;

    typedef /*typename*/ FunctionSpaceType::RangeFieldType RangeFieldType;


    // discrete function space
    typedef DiscreteFunctionSpace::Discontinuous::Orthonormal<FunctionSpaceType, GridPartType, polOrder> DGSpaceType;

    const DGSpaceType dgSpace(gridPart);


    //    // local evaluation
    //    typedef Dune::Functionals::Evaluation::Unary::Scale< FunctionSpaceType >
    //      ProductEvaluationType;

    //    ProductEvaluationType productEvaluation( "[1.0;1.0;1.0]", 0 );

    //    typedef Dune::Functionals::Evaluation::Binary::Elliptic< FunctionSpaceType >
    //      EllipticEvaluationType;

    //    EllipticEvaluationType ellipticEvaluation( "[1.0;1.0;1.0]", 0 );


    //    // operator and functional
    //    typedef DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
    //      LocalEllipticOperatorType;

    //    const LocalEllipticOperatorType localEllipticOperator( ellipticEvaluation );

    //    typedef DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
    //      LocalL2FunctionalType;

    //    const LocalL2FunctionalType localL2Functional( productEvaluation );


    //    // matrix, rhs and solution storage
    //    typedef Container::Matrix::Defaults< RangeFieldType, dimRange, dimRange >::BCRSMatrix
    //      MatrixFactory;

    //    typedef /*typename*/ MatrixFactory::AutoPtrType
    //      MatrixPtrType;

    //    MatrixPtrType A = MatrixFactory::create( discreteH1, discreteH1 );

    //    typedef Container::Vector::Defaults< RangeFieldType, dimRange >::BlockVector
    //      VectorFactory;

    //    typedef /*typename*/ VectorFactory::AutoPtrType
    //      VectorPtrType;

    //    VectorPtrType F = VectorFactory::create( discreteH1 );
    //    *F = 0.0;

    //    VectorPtrType u = VectorFactory::create( discreteH1 );
    //    *u = 0.0;


    //    // assembler
    //    typedef Assembler::Local::Codim0::Matrix< LocalEllipticOperatorType >
    //      LocalMatrixAssemblerType;

    //    const LocalMatrixAssemblerType localMatrixAssembler( localEllipticOperator );

    //    typedef Assembler::Local::Codim0::Vector< LocalL2FunctionalType >
    //      LocalVectorAssemblerType;

    //    const LocalVectorAssemblerType localVectorAssembler( localL2Functional );

    //    typedef Assembler::System::Unconstrained< DiscreteH1Type, DiscreteH1Type >
    //      SystemAssemblerType;

    //    const SystemAssemblerType systemAssembler( discreteH1, discreteH1 );

    //    systemAssembler.assembleSystem( localMatrixAssembler, *A,
    //                                    localVectorAssembler, *F );


    //    // preconditioner and solver
    //    typedef /*typename*/ MatrixFactory::ContainerType
    //      MatrixContainerType;

    //    typedef /*typename*/ VectorFactory::ContainerType
    //      VectorContainerType;

    //    typedef Dune::MatrixAdapter< MatrixContainerType, VectorContainerType, VectorContainerType >
    //      MatrixAdapterType;

    //    MatrixAdapterType matrix( *A );

    //    typedef Dune::SeqILU0< MatrixContainerType, VectorContainerType, VectorContainerType, 1 >
    //      PreconditionerType;

    //    PreconditionerType preconditioner( *A, 1.0 );

    //    typedef Dune::CGSolver< VectorContainerType >
    //      SolverType;

    //    SolverType solver( matrix, preconditioner, 1e-4, 100, 2 );

    //    Dune::InverseOperatorResult result;

    //    // u_0 = A^(-1) F
    //    solver.apply( *u, *F, result );


    //    // postprocessing
    //    typedef /*typename*/ Dune::Functionals::DiscreteFunction::Continuous::BlockVector< DiscreteH1Type >
    //      DiscreteFunctionType;

    //    DiscreteFunctionType solution( discreteH1, *u, "solution" );
    //    Dune::FemTools::DiscreteFunction::IO::writeToVTK( solution, "solution" );


    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
