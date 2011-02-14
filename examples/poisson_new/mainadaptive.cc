/**************************************************************************
**       Title: mainadaptive.cc
** Description: File demonstrating a simple numerics problem on arbitrary
**              grids: poisson-problem with known solution is given
**              and compared with numerical solution (EOC)
**              Dune grid parser is used.
**
**              For changing grid-types, compile with
**
**              make clean
**
**              and one of the following
**           
**              make                         or
**              make GRIDTYPE=YASPGRID       (default)
**                    -> compiles and works correctly with restrictions:
**                       - Adaptive scheme does not work, because YASPGRID
**                         doesn't implement conforming refinements
**                       - Parallel scheme only works at refinement level 0,
**                         because YASPGRID does not implement load balancing 
**              make GRIDTYPE=ALBERTAGRID
**                    -> compiles and works correctly
**                       - no parallelization implemented
**              make GRIDTYPE=SGRID
**                    -> compiles and works correctly
**                       - no parallelization implemented
**                       - no adaptivity implemented
**              make GRIDTYPE=ALUGRID_SIMPLEX
**                    -> compiles and works correctly
**                       - Adaptive scheme does not work, because
**                         ALUGRID_SIMPLEX doesn't implement conforming
**                         refinements. Use ALUGRID_CONFORM instead.
**
**
**              Similarly, the polynomial order can be specified by
**
**              make POLORDER=2
**
**************************************************************************/
#include <config.h>

// definition of algorithm
#include "algorithmadaptive.hh"

// definition of problem
#include "problemdata.hh"

// mpi manager calling MPI_Init (addon to mpihelper from DUNE)
#include <dune/fem/misc/mpimanager.hh>

#ifndef POLORDER 
const int polynomialOrder = 1;
#else 
const int polynomialOrder = POLORDER;
#endif

#include <dune/fem-howto/baseadaptive.hh>

typedef Dune::ProblemInterface< Dune::FunctionSpace< double, double, Dune::GridSelector::dimworld, 1 > > ProblemType;

// Main Program
// ------------                                                     /*@LST0S@*/
int main( int argc, char **argv )
{
  typedef Dune::GridSelector::GridType GridType;

  // initialize MPI 
  Dune::MPIManager :: initialize ( argc, argv );
  const int rank = Dune::MPIManager :: rank ();
  
  try
  {                                                                 /*@LST0E@*/
    // append parameters from the comand line 
    Dune::Parameter::append( argc, argv );

    // append parameters from the parameter file  
    Dune::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

    // generate GridPointer holding grid instance                   /*@LST0S@*/
    Dune::GridPtr< GridType > gridptr = initialize< GridType >(std::string("Poisson problem"));/*@\label{poi:gridinit}@*/

    // get grid reference 
    GridType& grid = *gridptr ;
      
    // create problem 
    ProblemType* problem = Dune::createProblem<GridType> ();
    assert( problem );

    // create stepper class 
    AdaptiveAlgorithm<GridType, polynomialOrder> algorithm(grid, *problem);  /*@\label{poi:stepper}@*/
    // compute solution
    computeAdaptive(algorithm);                                       /*@\label{poi:compute}@*//*@LST0E@*/

    // write parameter logfile 
    Dune::Parameter :: write("parameter.log");

    // remove problem 
    delete problem;

    return 0;
  }                                                                 /*@LST0S@*/
  catch( Dune::Exception &exception )                                     /*@\label{poi:catch0}@*/
  {
    if( rank == 0 )
      std :: cerr << exception << std :: endl;
    return 1;
  }                                                                 /*@\label{poi:catch1}@*//*@LST0E@*/
}
