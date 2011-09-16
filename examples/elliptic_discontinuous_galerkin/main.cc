/**
  \file   main.cc
  \brief  Main file for the elliptic discontinuous galerkin example.
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
//#include <dune/common/fmatrix.hh>

// dune-grid includes
#include <dune/grid/utility/gridtype.hh>

// dune-istl includes
#include <dune/istl/solvers.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/misc/mpimanager.hh>

// reenable warnings
#include <dune/helper-tools/header/enablewarnings.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/string.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed-discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim1/integral.hh>
#include <dune/detailed-discretizations/container/factory.hh>
//#include <dune/detailed-discretizations/assembler/local/codim1/matrix.hh>

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

    const std::string dgfFilename = "../macrogrids/unitcube"
      + Dune::HelperTools::Common::String::toString( GRIDDIM ) + ".dgf";

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
      DGSpaceType;

    const DGSpaceType dgSpace( gridPart );


    // local evaluation
//    typedef Dune::Functionals::Evaluation::Local::Unary::Scale< FunctionSpaceType >
//      ProductEvaluationType;

//    const ProductEvaluationType productEvaluation( "[1.0;1.0;1.0]", 0 );

//    typedef Dune::Functionals::Evaluation::Loca::Binary::Elliptic< FunctionSpaceType >
//      EllipticEvaluationType;

//    const EllipticEvaluationType ellipticEvaluation( "[1.0;1.0;1.0]", 0 );

    typedef Dune::DetailedDiscretizations::Evaluation::Local::Quaternary::IPDGFluxes::Inner< FunctionSpaceType >
      InnerFacesIPDGEvaluationType;

    const InnerFacesIPDGEvaluationType innerFacesIPDGEvaluation( "[1.0;1.0;1.0]" );

    typedef Dune::DetailedDiscretizations::Evaluation::Local::Quaternary::IPDGFluxes::Dirichlet< FunctionSpaceType >
      DirichletFacesIPDGEvaluationType;

    const DirichletFacesIPDGEvaluationType dirichletFacesIPDGEvaluation( "[1.0;1.0;1.0]" );


    // operator and functional
//    typedef DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
//      LocalEllipticOperatorType;

//    const LocalEllipticOperatorType localEllipticOperator( ellipticEvaluation );

//    typedef DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
//      LocalL2FunctionalType;

//    const LocalL2FunctionalType localL2Functional( productEvaluation );

    typedef DiscreteOperator::Local::Codim1::Integral< InnerFacesIPDGEvaluationType >
      LocalIPDGInnerFacesOperatorType;

    const LocalIPDGInnerFacesOperatorType localIPDGInnerFacesOperator( innerFacesIPDGEvaluation );


    // matrix, rhs and solution storage
    typedef Container::Matrix::Defaults< RangeFieldType, dimRange, dimRange >::BCRSMatrix
      MatrixFactory;

    typedef /*typename*/ MatrixFactory::AutoPtrType
      MatrixPtrType;

    MatrixPtrType A = MatrixFactory::create( dgSpace, dgSpace );

    typedef Container::Vector::Defaults< RangeFieldType, dimRange >::BlockVector
      VectorFactory;

    typedef /*typename*/ VectorFactory::AutoPtrType
      VectorPtrType;

    VectorPtrType F = VectorFactory::create( dgSpace );
    *F = 0.0;

    VectorPtrType u = VectorFactory::create( dgSpace );
    *u = 0.0;


    // assembler
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

