/**
  \file   main.cc
  \brief  Main file for the elliptic discontinuous galerkin example.
  **/

// disable warnings about problems in dune headers
#include <dune/helper-tools/header/disablewarnings.hh>

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
#include <dune/fem/space/common/functionspace.hh>
//#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/misc/mpimanager.hh>

// reenable warnings
#include <dune/helper-tools/header/enablewarnings.hh>

// dune-helper-tools includes
#include <dune/helper-tools/common/string.hh>
#include <dune/helper-tools/discretefunction/io.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/discretefunctionspace/discontinuous/orthonormal.hh>
#include <dune/detailed-discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed-discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed-discretizations/evaluation/local/quaternary/ipdgfluxes.hh>
#include <dune/detailed-discretizations/evaluation/local/unary/ipdgfluxes.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim1/innerintegral.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim1/boundaryintegral.hh>
#include <dune/detailed-discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed-discretizations/discretefunctional/local/codim1/integral.hh>
#include <dune/detailed-discretizations/discretefunction/discontinuous.hh>
#include <dune/detailed-discretizations/container/factory.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed-discretizations/assembler/local/codim1/matrix.hh>
#include <dune/detailed-discretizations/assembler/local/combined/matrix.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed-discretizations/assembler/local/codim1/vector.hh>
#include <dune/detailed-discretizations/assembler/local/combined/vector.hh>
#include <dune/detailed-discretizations/assembler/system/unconstrained.hh>

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

    typedef typename GridType::LeafGridView
      GridViewType;

    const std::string dgfFilename = "../macrogrids/unitcube"
      + Dune::HelperTools::Common::String::toString( GRIDDIM ) + ".dgf";

    Dune::GridPtr< GridType > gridPtr( dgfFilename );

    const GridType& grid( *gridPtr );

    const GridViewType& gridView = grid.leafView();

    static const unsigned int dimDomain = GridType::dimension;

    static const unsigned int dimRange = 1;

    // function space
    typedef Dune::FunctionSpace< double, double, dimDomain, dimRange >
      FunctionSpaceType;

    typedef typename FunctionSpaceType::RangeFieldType
      RangeFieldType;


    // discrete function space
    typedef DiscreteFunctionSpace::Discontinuous::Orthonormal< FunctionSpaceType, GridViewType, polOrder >
      DGSpaceType;

    const DGSpaceType dgSpace( gridView );


    // local evaluation
    typedef Evaluation::Local::Unary::Scale< FunctionSpaceType >
      ProductEvaluationType;

    const ProductEvaluationType productEvaluation( "[1.0;1.0;1.0]", 0 );

    typedef Evaluation::Local::Binary::Elliptic< FunctionSpaceType >
      EllipticEvaluationType;

    const EllipticEvaluationType ellipticEvaluation( "[1.0;1.0;1.0]", 0 );

    typedef Evaluation::Local::Quaternary::IPDGFluxes::Inner< FunctionSpaceType >
      InnerFacesIPDGEvaluationType;

    const InnerFacesIPDGEvaluationType innerFacesIPDGEvaluation( "[1.0;1.0;1.0]" );

    typedef Evaluation::Local::Quaternary::IPDGFluxes::Dirichlet< FunctionSpaceType >
      DirichletFacesIPDGEvaluationType;

    const DirichletFacesIPDGEvaluationType dirichletFacesIPDGEvaluation( "[1.0;1.0;1.0]", 0 );

    typedef Evaluation::Local::Unary::IPDGFluxes< FunctionSpaceType >
      IPDGEvaluationType;

    const IPDGEvaluationType ipdgEvaluation(  "[1.0;1.0;1.0]",
                                              "[0.0;0.0;0.0]",
                                              0 );


    // operator and functional
    typedef DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType >
      LocalEllipticOperatorType;

    const LocalEllipticOperatorType localEllipticOperator( ellipticEvaluation );

    typedef DiscreteOperator::Local::Codim1::InnerIntegral< InnerFacesIPDGEvaluationType >
      LocalIPDGInnerFacesOperatorType;

    const LocalIPDGInnerFacesOperatorType localIPDGInnerFacesOperator( innerFacesIPDGEvaluation );

    typedef DiscreteOperator::Local::Codim1::BoundaryIntegral< DirichletFacesIPDGEvaluationType >
      LocalIPDGDirichletFacesOperatorType;

    const LocalIPDGDirichletFacesOperatorType localIPDGDirichletFacesOperator( dirichletFacesIPDGEvaluation );

    typedef DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType >
      LocalL2FunctionalType;

    const LocalL2FunctionalType localL2Functional( productEvaluation );

    typedef DiscreteFunctional::Local::Codim1::Integral< IPDGEvaluationType >
      LocalIPDGFunctionalType;

    const LocalIPDGFunctionalType localIPDGFunctional( ipdgEvaluation );


    // matrix, rhs and solution storage
    typedef Container::Matrix::Defaults< RangeFieldType, dimRange, dimRange >::BCRSMatrix
      MatrixFactory;

    typedef typename MatrixFactory::AutoPtrType
      MatrixPtrType;

    MatrixPtrType A = MatrixFactory::create( dgSpace, dgSpace );

    typedef Container::Vector::Defaults< RangeFieldType, dimRange >::BlockVector
      VectorFactory;

    typedef typename VectorFactory::AutoPtrType
      VectorPtrType;

    VectorPtrType F = VectorFactory::create( dgSpace );
    *F = 0.0;

    VectorPtrType u = VectorFactory::create( dgSpace );
    *u = 0.0;


    // local matrix assembler
    typedef Assembler::Local::Codim0::Matrix< LocalEllipticOperatorType >
      LocalCodim0MatrixAssemblerType;

    const LocalCodim0MatrixAssemblerType localCodim0MatrixAssembler( localEllipticOperator );

    typedef Assembler::Local::Codim1::Matrix< LocalIPDGInnerFacesOperatorType, LocalIPDGDirichletFacesOperatorType >
      LocalCodim1MatrixAssemblerType;

    const LocalCodim1MatrixAssemblerType localCodim1MatrixAssembler(  localIPDGInnerFacesOperator,
                                                                      localIPDGDirichletFacesOperator );

    typedef Assembler::Local::Combined::Matrix< LocalCodim0MatrixAssemblerType, LocalCodim1MatrixAssemblerType >
      LocalCombinedMatrixAssemblerType;

    const LocalCombinedMatrixAssemblerType localCombinedMatrixAssembler(  localCodim0MatrixAssembler,
                                                                          localCodim1MatrixAssembler );


    // local vector assembler
    typedef Assembler::Local::Codim0::Vector< LocalL2FunctionalType >
      LocalCodim0VectorAssemblerType;

    const LocalCodim0VectorAssemblerType localCodim0VectorAssembler( localL2Functional );

    typedef Assembler::Local::Codim1::Vector< LocalIPDGFunctionalType >
      LocalCodim1VectorAssemblerType;

    const LocalCodim1VectorAssemblerType localCodim1VectorAssembler( localIPDGFunctional );

    typedef Assembler::Local::Combined::Vector< LocalCodim0VectorAssemblerType, LocalCodim1VectorAssemblerType >
      LocalCombinedVectorAssemblerType;

    const LocalCombinedVectorAssemblerType localCombinedVectorAssembler(  localCodim0VectorAssembler,
                                                                          localCodim1VectorAssembler );


    // system assembler
    typedef Assembler::System::Unconstrained< DGSpaceType, DGSpaceType >
      SystemAssemblerType;

    const SystemAssemblerType systemAssembler( dgSpace, dgSpace );

    systemAssembler.assembleSystem( localCombinedMatrixAssembler, *A,
                                    localCombinedVectorAssembler, *F );


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

    // u = A^(-1) F
    solver.apply( *u, *F, result );


    // postprocessing
    typedef typename DiscreteFunction::Discontinuous::BlockVector< DGSpaceType >
      DiscreteFunctionType;

    DiscreteFunctionType solution( dgSpace, *u, "solution" );
    Dune::HelperTools::DiscreteFunction::IO::writeToVTK( solution, "solution" );


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

