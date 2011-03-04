#ifndef POISSON_SERIAL_HH
#define POISSON_SERIAL_HH

//***********************************************************************
/*! Poisson problem: 

  This is an example solving the poisson problem
  \f{displaymath}
  \begin{array}{rcll}
  -\triangle u &=& f & \quad \mbox{in}\ \Omega\\
             u &=& 0 & \quad \mbox{on}\ \partial\Omega
  \end{array}
  \f}
  with the finite element method using Lagrangian elements. The polynomial
  order is given by POLORDER.
  
  In this example, $\Omega = ]0,1[^{dimworld}$ and
  \f[
  f( x, y, z ) = 4 dimworld \pi^2 u(x,y,z) 
  \f]

  The exact solution to the poisson problem is then given by
  \f[
  u( x, y, z ) = \prod_{i=1}^{dimworld} sin( 2 \pi x_i).
  \f]
*/
//***********************************************************************

// this define enables the use of higher (>2) order Lagrange basis functions
#define USE_TWISTFREE_MAPPER

// System Includes
// ---------------
#include <iostream>
#include <sstream>


// DUNE Core Includes
// ------------------
#include <dune/common/version.hh>

// if this value is not defined, then we have version 1.1.1
#ifndef DUNE_VERSION_HH
#define OLD_DUNE_GRID_VERSION
#endif

#include <dune/common/stdstreams.hh>
#include <dune/common/timer.hh>

// DUNE-FEM includes
// -----------------

// grid parts 
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// adaptation classes 
#include <dune/fem/space/common/adaptmanager.hh>
// lagrange space 
#include <dune/fem/space/lagrangespace.hh>

// discrete functions 
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/vectorfunction.hh>

#include <dune/fem/operator/discreteoperatorimp.hh>

// matrix implementations 
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/operator/matrix/ontheflymatrix.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>

// linear solvers 
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/solver/istlsolver.hh>

// l2 norm and h1 norm 
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// parameter reading facilities 
#include <dune/fem/io/parameter.hh>

// header for data output for grape and VTK 
#include <dune/fem/io/file/datawriter.hh>

#include <dune/fem/misc/femeoc.hh>

// Local Includes
// --------------

// matrix assembler 
#include "laplaceoperator.hh"
#include "functional.hh"
#include "dirichletconstraints.hh"
#include "neumannconstraints.hh"


// Algorithm
// ---------

template <class GridImp, int polynomialOrder,
          class GridPartImp = Dune::AdaptiveLeafGridPart<GridImp,Dune::InteriorBorder_Partition> >
          //class GridPartImp = Dune::LeafGridPart<GridImp> >
class Algorithm                                                       /*@LST0S@*/
{
public:                                                               /*@LST0E@*/
  typedef GridImp HGridType;

  /** choose the grid partition (and hence the index set) to use
   *
   *  \note Not all index sets are continuous. The LeafIndexSet for AlbertaGrid,
   *        for example, is not. If you want to use OEM solvers, the index set
   *        must be continuous. In such a case use AdaptiveLeafGridPart.
   */
  //---- GridParts -----------------------------------------------------------
  typedef GridPartImp  GridPartType;
                                                                    
  //---- FunctionSpace -------------------------------------------------------
  //! define the function space, \f[ \R^n \rightarrow \R \f]
  // see dune/common/functionspace.hh
  typedef Dune::FunctionSpace< double, double, HGridType::dimensionworld, 1 > FunctionSpaceType;

  //! to be revised 
  typedef Dune::ProblemInterface< FunctionSpaceType > ProblemType;

  // The data functions (as defined in problem.hh)
  //---- Right Hand Side, Exact Solution, and Stiffness tensor ---------------
  typedef typename ProblemType ::  ExactSolutionType   ExactSolutionType;

  //---- Adapter for exact solution ------------------------------------------
  typedef Dune::DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
    GridExactSolutionType;

  //---- DiscreteFunctionSpace -----------------------------------------------
  //! define the discrete function space our unkown belongs to
  typedef Dune::LagrangeDiscreteFunctionSpace
    < FunctionSpaceType, GridPartType, polynomialOrder, Dune::CachingStorage >
    DiscreteSpaceType;

  //---- DiscreteFunction ----------------------------------------------------
  //---- good choice for adaptive simulations using OEM solver ---------------
  //! define the type of discrete function we are using
  typedef Dune::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
  //---- other possible choices, use BlockVectorDiscreteFunction for ISTL ----
  //typedef Dune::BlockVectorDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
  //typedef Dune::ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteSpaceType, DynamicVector< double > > > DiscreteFunctionType;

  //---- MatrixObjects -------------------------------------------------------
  //---- good choice for build in solvers ------------------------------------
  //! define the type of the system matrix object
  typedef Dune::SparseRowMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > 
      MatrixObjectTraits;

  //---- other choices, ISTLMatrixTraits for BCRSMatrix from DUNE-ISTL -------
  //typedef Dune::ISTLMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > MatrixObjectTraits;
  //typedef Dune::OnTheFlyMatrixTraits < DiscreteSpaceType, DiscreteSpaceType > MatrixObjectTraits;
      
  //! define the discrete laplace operator, see ./laplaceoperator.hh 
  typedef Dune::LaplaceOperator< DiscreteFunctionType, MatrixObjectTraits >
    LaplaceOperatorType;                                            

  //---- InverseOperator ----------------------------------------------------
  //---- good choice for build in CG solver ---------------------------------
  //! define the inverse operator we are using to solve the system 
  typedef Dune::CGInverseOp< DiscreteFunctionType, LaplaceOperatorType >
     InverseOperatorType;
    //---- other choices ------------------------------------------------------ 
    //typedef Dune::ISTLBICGSTABOp< DiscreteFunctionType, LaplaceOperatorType >
    //typedef Dune::ISTLCGOp< DiscreteFunctionType, LaplaceOperatorType >
    //typedef Dune::OEMCGOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef Dune::OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef Dune::OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType>
    //typedef Dune::OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType>
    //    InverseOperatorType;

public:
  //! constructor                                                  /*@LST0S@*/
  Algorithm(HGridType & grid, const ProblemType& problem)
  : grid_( grid ),
    gridPart_( grid_ ),
    space_( gridPart_ ),
    problem_( problem )
  {
    // add entries to eoc calculation 
    std::vector<std::string> femEocHeaders;
    // we want to calculate L2 and H1 error (and EOC)
    femEocHeaders.push_back("$L^2$-error");
    femEocHeaders.push_back("$H^1$-error");                        /*@LST0E@*/

    // get eoc id 
    eocId_ = Dune::FemEoc::addEntry(femEocHeaders);

    // distribute grid (only for parallel runs, in serial runs nothing is done here)
    grid.loadBalance();
  }                                                                /*@LST0S@*/

  //! setup and solve the linear system
  void operator()(DiscreteFunctionType & solution)
  {                                                                /*@LST0E@*/
    // in verbose mode some info 
    if( Dune::Parameter :: verbose () ) 
    {
      std::cout << std :: endl << "Solving for " << space_.size()
        << " unkowns and polynomial order "
        << DiscreteSpaceType :: polynomialOrder << "." 
        << std :: endl << std :: endl;
    }

    // initialize solution with zero                                 /*@LST0S@*/
    solution.clear();

    // create laplace assembler (is assembled on call of systemMatrix by solver)
    LaplaceOperatorType laplace( space_ , problem_ );               /*@\label{poi:laplaceDef}@*/

    // functional describing the right hand side 
    typedef Dune::IntegralFunctional< ProblemType, DiscreteFunctionType > FunctionalType ;
    FunctionalType rhsFunctional( problem_ );

    Dune::DirichletConstraints< DiscreteSpaceType, ProblemType > 
      constraints( space_, problem_ );                              /*@\label{poi:rhsInit1}@*/

    Dune::NeumannConstraints< DiscreteSpaceType, ProblemType >
      neumann( space_, problem_ );

    // solve the linear system  @\ref{
    solve( laplace, rhsFunctional, solution , constraints, neumann );  /*@\label{poi:solve}@*/ 
  }

  //! finalize computation by calculating errors and EOCs 
  void finalize(DiscreteFunctionType & solution)
  {
    // create exact solution 
    ExactSolutionType uexact( problem_ ); 

    // create grid function adapter 
    GridExactSolutionType ugrid( "exact solution", uexact, gridPart_,
                                 DiscreteSpaceType :: polynomialOrder + 1 );

    // create L2 - Norm 
    Dune::L2Norm< GridPartType > l2norm( gridPart_ );
    // calculate L2 - Norm 
    const double l2error = l2norm.distance( ugrid, solution );      

    // create H1 - Norm 
    Dune::H1Norm< GridPartType > h1norm( gridPart_ );
    // calculate H1 - Norm 
    const double h1error = h1norm.distance( ugrid, solution );

    // store values 
    std::vector<double> errors;
    errors.push_back( l2error );
    errors.push_back( h1error );

    // submit error to the FEM EOC calculator 
    Dune::FemEoc :: setErrors(eocId_, errors);
  }                                                                 

  //! return reference to discrete space 
  DiscreteSpaceType & space() { return space_; }                    /*@LST0E@*/

private:                                                            /*@LST0S@*/
  //! solve the resulting linear system
  template <class FunctionalType, 
            class Constraints,
            class Neumann >
  void solve ( LaplaceOperatorType &laplace,
               const FunctionalType& functional,
               DiscreteFunctionType &solution,
               const Constraints& constraints,
               const Neumann& neumann )
  {                                                                /*@LST0E@*/
    // create storage for rhs (this could be any vector type)
    DiscreteFunctionType rhs( "rhs", space_ );

    // assemble right hand side functional into vector 
    Dune::AssembledFunctional< FunctionalType > rhsFunctional ( space_, functional );
    rhsFunctional.assemble( rhs );

    // apply constraints (here Dirichlet boundary) 
    bool hasDirichletBoundary = 
      constraints.apply( laplace.systemMatrix(), rhs, solution );

    if( ! hasDirichletBoundary ) 
      neumann.apply( laplace.systemMatrix(), rhs, solution );

    // create timer (also stops time) 
    Dune::Timer timer;

    // solve the linear system (with CG)
    const double reduction = 1e-6;
    double solverEps = 1e-8 ;
    solverEps = Dune::Parameter :: getValue( "femhowto.solvereps", solverEps );

    // get verbosity information from Parameter class 
    const bool verbose = Dune::Parameter :: verbose();

    // create inverse operator (note reduction is not used here)    /*@LST0S@*/
    InverseOperatorType cg( laplace, reduction, solverEps );/*@\label{poi:cgInit}@*/

    // solve the system 
    cg( rhs, solution );                                            /*@\label{poi:cgExec}@*//*@LST0E@*/

    // get elapsed time 
    const double solveTime = timer.elapsed();

    // output in only in verbose mode (parameter fem.verbose)
    if( verbose ) 
      std :: cout << "Time needed by solver: " << solveTime << "s" << std :: endl;
  }                                                                /*@LST0S@*/

protected:
  HGridType& grid_;            // reference to grid, i.e. the hierarchical grid 
  GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid 
  DiscreteSpaceType space_;    // the discrete function space
  const ProblemType& problem_; // the problem data
  int eocId_;                  // id for FemEOC
};                                                                /*@LST0E@*/

#endif // POISSONSERIAL_HH
