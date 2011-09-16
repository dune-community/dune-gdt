/**
  \file   main.cc
  \brief  Main file for the elliptic finite element example.
  **/

// disable warnings about problems in dune headers
#include <dune/helper-tools/header/disablewarnings.hh>

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
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/common/functionspace.hh>

// reenable warnings
#include <dune/helper-tools/header/enablewarnings.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed-discretizations/discretefunctionspace/subspace/linear.hh>
#include <dune/detailed-discretizations/discretefunctionspace/subspace/affine.hh>
#include <dune/detailed-discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed-discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed-discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed-discretizations/container/factory.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed-discretizations/assembler/system/affine.hh>
#include <dune/detailed-discretizations/discretefunction/continuous.hh>
#include <dune/detailed-discretizations/discretefunction/femadapter.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/string.hh>
#include <dune/helper-tools/discretefunction/io.hh>

// only ever do this in a .cc!
using namespace Dune::DetailedDiscretizations;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif


int main( int argc, char** argv )
{
  try{

    // MPI manager
    Dune::MPIManager::initialize( argc, argv );


    // grid
    typedef Dune::GridSelector::GridType
      GridType;

    typedef Dune::LeafGridPart< GridType >
      GridPartType;

    const std::string dgfFilename = "../macrogrids/unitcube" +
      Dune::HelperTools::Common::String::toString( GRIDDIM ) + ".dgf";

    Dune::GridPtr< GridType > gridPtr( dgfFilename );

    GridPartType gridPart( *gridPtr );

    static const unsigned int dimDomain = GridType::dimension;

    static const unsigned int dimRange = 1;


    // function space
    typedef Dune::FunctionSpace< double, double, dimDomain, dimRange >
      FunctionSpaceType;

    typedef /*typename*/ FunctionSpaceType::RangeFieldType
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

    const DiscreteH1GType discreteH1G( discreteH10, "[1.0+0.01x+0.01y;y;z]" );


    // local evaluation
    typedef Evaluation::Local::Unary::Scale< FunctionSpaceType >
      ProductEvaluationType;

    ProductEvaluationType productEvaluation( "[1.0;1.0;1.0]", 0 );

    typedef Evaluation::Local::Binary::Elliptic< FunctionSpaceType >
      EllipticEvaluationType;

    EllipticEvaluationType ellipticEvaluation( "[1.0;1.0;1.0]", 0 );


    // operator and functional
    typedef DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
      LocalEllipticOperatorType;

    const LocalEllipticOperatorType localEllipticOperator( ellipticEvaluation );

    typedef DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
      LocalL2FunctionalType;

    const LocalL2FunctionalType localL2Functional( productEvaluation );


    // matrix, rhs and solution storage
    typedef Container::Matrix::Defaults< RangeFieldType, dimRange, dimRange >::BCRSMatrix
      MatrixFactory;

    typedef /*typename*/ MatrixFactory::AutoPtrType
      MatrixPtrType;

    MatrixPtrType A = MatrixFactory::create( discreteH1, discreteH1 );

    typedef Container::Vector::Defaults< RangeFieldType, dimRange >::BlockVector
      VectorFactory;

    typedef /*typename*/ VectorFactory::AutoPtrType
      VectorPtrType;

    VectorPtrType F = VectorFactory::create( discreteH1 );
    *F = 0.0;

    VectorPtrType G = VectorFactory::create( discreteH1 );
    *G = 0.0;

    VectorPtrType u = VectorFactory::create( discreteH1 );
    *u = 0.0;


    // assembler
    typedef Assembler::Local::Codim0::Matrix< LocalEllipticOperatorType >
      LocalMatrixAssemblerType;

    const LocalMatrixAssemblerType localMatrixAssembler( localEllipticOperator );

    typedef Assembler::Local::Codim0::Vector< LocalL2FunctionalType >
      LocalVectorAssemblerType;

    const LocalVectorAssemblerType localVectorAssembler( localL2Functional );

    typedef Assembler::System::Affine< DiscreteH1GType, DiscreteH10Type >
      SystemAssemblerType;

    const SystemAssemblerType systemAssembler( discreteH1G, discreteH10 );

    systemAssembler.assembleSystem( localMatrixAssembler, *A,
                                    localVectorAssembler, *F,
                                    *G );


    // preconditioner and solver
    typedef /*typename*/ MatrixFactory::ContainerType
      MatrixContainerType;

    typedef /*typename*/ VectorFactory::ContainerType
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
    *F -= *G;
    solver.apply( *u, *F, result );

    // u = u0 + g
    *u += discreteH1G.affineShift().storage();


    // postprocessing
    typedef /*typename*/ DiscreteFunction::Continuous::BlockVector< DiscreteH1Type >
      DiscreteFunctionType;

    DiscreteFunctionType solution( discreteH1, *u, "solution" );

    typedef /*typename*/ DiscreteFunction::FemAdapter< DiscreteFunctionType >
      DiscreteFunctionFemAdapterType;

    DiscreteFunctionFemAdapterType solutionFemAdapter( solution );

    Dune::HelperTools::DiscreteFunction::IO::writeToVTK( solutionFemAdapter, "solution" );


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
