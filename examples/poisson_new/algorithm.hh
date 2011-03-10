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
#include <string>

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

#include <dune/grid/common/gridenums.hh>

// Local Includes
// --------------

// matrix assembler
#include "laplaceoperator.hh"
#include "functional.hh"
#include "dirichletconstraints.hh"
#include "neumannconstraints.hh"

//new constraints implementation
//#include "constraints.hh"


// Algorithm
//------------



/**
 * @class LinearSubspace
 *
 * @brief A class representing a linear subspace.
 *
 * @tparam DiscreteFunctionSpaceImp The discrete function space where the sub space is constructed from.
 * @tparam ConstraintsImp The constraints needed to construct a sub space. 
 * Use for example Dune::Functionals::DirichletConstraints or Dune::Functionals::NeumannConstraints.
 *
 * @todo improve comments
 *
 * @attention class is not finished yet
 */
template< class DiscreteFunctionSpaceImp, class ConstraintsImp >
class LinearSubspace
  : public DiscreteFunctionSpaceImp
{
public:
  typedef ConstraintsImp
    ConstraintsType;
  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::Traits
    Traits;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType
    GridType;
  typedef typename DiscreteFunctionSpaceType::IndexSetType
    IndexSetType;
  typedef typename DiscreteFunctionSpaceType::IteratorType
    IteratorType;
  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;
  typedef typename DiscreteFunctionSpaceType::DofManagerType
    DofManagerType;
  typedef typename DiscreteFunctionSpaceType::CommunicationManagerType
    CommunicationManagerType;
  typedef typename DiscreteFunctionSpaceType::MapperType
    MapperType;
  typedef typename DiscreteFunctionSpaceType::BlockMapperType
    BlockMapperType;
  enum {localBlockSize = DiscreteFunctionSpaceType::localBlockSize};
  //enum {polynomialOrder = DiscreteFunctionType::polynomialOrder};

private:
  typedef DiscreteFunctionSpaceType
    BaseType;

public:
  /**
   * @brief constructor
   *
   * @param space The linear subspace.
   * @param constraints The constraints.
   * @param commInterface interface type
   * @param commDirection communication direction for parallel runs
   */
  LinearSubspace( DiscreteFunctionSpaceType& space,
                  const ConstraintsType& constraints,
                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication )
    : BaseType( space.gridPart(), commInterface, commDirection ), constraints_( constraints )
  {
  }

  /**
   * @brief constructor
   *
   * @param gridPart Grid part.
   * @param constraints The constraints.
   * @param commInterface interface type
   * @param commDirection communication direction for parallel runs
   */
  LinearSubspace( GridPartType& gridPart,
                  const ConstraintsType& constraints,
                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication )
    : BaseType( gridPart, commInterface, commDirection ), constraints_( constraints )
  {
  }
  /**
   * @brief get constraints
   */
  const ConstraintsType getConstraints() const
  {
    return constraints_;
  }

private:
  const ConstraintsType& constraints_;

};


/**
 * @class LinearSubspaceDiscreteFunction
 *
 * @brief A class representing a discrete function from a linear subspace.
 *
 * @tparam LinearSubspaceImp A linear subspace.
 * @tparam DiscreteFunctionImp A discrete function where the linear sub space is constructed from.
 *
 * @todo improve comments
 *
 * @attention class is not finished yet
 */
template< class LinearSubspaceImp, class DiscreteFunctionImp >
class LinearSubspaceDiscreteFunction
  : public DiscreteFunctionImp
{
public:
  typedef DiscreteFunctionImp
    DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef LinearSubspaceImp
    LinearSubspaceType;

  //TODO: assert: LinearSubspaceType::DiscreteFunctionSpaceType and
  //              DiscreteFunctionType::DiscreteFunctionSpaceType should fit

  typedef typename DiscreteFunctionType::DomainFieldType
    DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType
    RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType
    DomainType;
  typedef typename DiscreteFunctionType::RangeType
    RangeType;

  typedef typename DiscreteFunctionType::JacobianRangeType
    JacobianRangeType;
  typedef typename DiscreteFunctionType::HessianRangeType
    HessianRangeType;

  typedef typename DiscreteFunctionType::GridType
    GridType;
  typedef typename DiscreteFunctionType::LocalFunctionType
    LocalFunctionType;
  typedef typename DiscreteFunctionType::DofIteratorType
    DofIteratorType;
  typedef typename DiscreteFunctionType::ConstDofIteratorType
    ConstDofIteratorType;
  typedef typename DiscreteFunctionType::MappingType
    MappingType;
  typedef typename DiscreteFunctionType::FunctionSpaceType
    FunctionSpaceType;
  typedef typename DiscreteFunctionType::FunctionType
    FunctionType;

private:
  typedef DiscreteFunctionImp
    BaseType;

public:

  /**
   * @brief constructor
   *
   * The constructor automatically applys the constraints
   *
   * @param discFunc Discrete function without any constraints living 
   * on the discrete function space the LinearSubspaceType was constructed from.
   * @param linDiscFuncSpace Discrete function space with constraints, i.e. a linear subspace.
   */
  LinearSubspaceDiscreteFunction( DiscreteFunctionType& discFunc,
                                  LinearSubspaceType& linDiscFuncSpace )
    : discFunc_( discFunc ),
      linDiscFuncSpace_( linDiscFuncSpace )
  {
     //TODO uncomment and correct it
     //ConstraintsType constraints = linDiscFuncSpace.getConstraints();
     //constraints.apply();
  }

  /**
   * @brief Returns the discrete function.
   */
  DiscreteFunctionType& discreteFunction()
  {
    return discFunc_;
  }

  /**
   * @brief Returns the discrete function space (here: linear subspace).
   */

  LinearSubspaceType& linearSubspace()
  {
    return linDiscFuncSpace_;
  }

private:
  DiscreteFunctionType& discFunc_;
  LinearSubspaceType& linDiscFuncSpace_;

};



/**
 * @class AffineSubspace
 *
 * @brief A class representing an affine subspace.
 *
 * We need this class for representing boundary data.
 *
 * @tparam LinearSubspaceImp The linear subspace.
 * @tparam DiscreteFunctionAffinePartImp The discrete function representing the affine part.
 *
 * @todo improve comments
 *
 * @attention class is not finished yet
 */
template< class LinearSubspaceImp, class DiscreteFunctionAffinePartImp >
class AffineSubspace
  : public LinearSubspaceImp
{
public:
  typedef DiscreteFunctionAffinePartImp
    DiscreteFunctionAffinePartType;
  typedef LinearSubspaceImp
    LinearSubspaceType;

  typedef typename LinearSubspaceType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename LinearSubspaceType::ConstraintsType
    ConstraintsType;

  typedef typename LinearSubspaceType::Traits
    Traits;

  typedef typename LinearSubspaceType::BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename LinearSubspaceType::GridPartType
    GridPartType;
  typedef typename LinearSubspaceType::GridType
    GridType;
  typedef typename LinearSubspaceType::IndexSetType
    IndexSetType;
  typedef typename LinearSubspaceType::IteratorType
    IteratorType;
  typedef typename GridPartType::IntersectionIteratorType
    IntersectionIteratorType;
  typedef typename LinearSubspaceType::EntityType
    EntityType;
  typedef typename LinearSubspaceType::DofManagerType
    DofManagerType;
  typedef typename LinearSubspaceType::CommunicationManagerType
    CommunicationManagerType;
  typedef typename LinearSubspaceType::MapperType
    MapperType;
  typedef typename LinearSubspaceType::BlockMapperType
    BlockMapperType;

private:
  typedef LinearSubspaceImp
    BaseType;
  typedef typename LinearSubspaceType::DomainType
    DomainType;
  typedef typename LinearSubspaceType::RangeType
    RangeType;
  typedef typename LinearSubspaceType::DomainFieldType
    DomainFieldType;
  typedef typename LinearSubspaceType::RangeFieldType
    RangeFieldType;

public:
  enum {localBlockSize = LinearSubspaceType::localBlockSize};
  
  /**
   * @brief constructor
   *
   * @param space The linear sub space.
   * @param df The affine part we want to "add" to the linear subspace
   * @param commInterface interface type
   * @param commDirection communication direction for parallel runs
   */
  AffineSubspace( LinearSubspaceImp& space, 
                  const DiscreteFunctionAffinePartType df,
                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication )
    : BaseType( space.gridPart(), space.getConstraints(), commInterface, commDirection ), df_( &df )
  {
  }

//  /**
//   * @brief constructor
//   *
//   * @param gridPart The grid part the function space is living on.
//   * @param df The affine part we want to "add" to the linear subspace
//   * @param commInterface interface type
//   * @param commDirection communication direction for parallel runs
//   */
//  AffineSubspace( GridPartType& gridPart, 
//                  const DiscreteFunctionAffinePartType df,
//                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
//                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication )
//    : BaseType( gridPart, commInterface, commDirection ), df_( &df )
//  {
//  }
//
//   /**
//   * @brief constructor
//   *
//   * @param space The linear subspace.
//   * @param df The affine part we want to "add" to the linear subspace
//   * @param commInterface interface type
//   * @param commDirection communication direction for parallel runs
//   */
//  AffineSubspace( LinearSubspaceType& space, 
//                  const DiscreteFunctionAffinePartType df,
//                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
//                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication )
//    : BaseType( space.gridPart(), commInterface, commDirection ), space_( space ), df_( &df )
//  {
//  }
//
  /**
   * @brief constructor
   *
   * @param gridPart The grid part the function space is living on.
   * @param commInterface interface type
   * @param commDirection communication direction for parallel runs
   */
  AffineSubspace( LinearSubspaceType& space, 
                  const Dune::InterfaceType commInterface = Dune::InteriorBorder_All_Interface,
                  const Dune::CommunicationDirection commDirection = Dune::ForwardCommunication)
    : BaseType( space.gridPart(), space.getConstraints(), commInterface, commDirection ), space_( space ), df_( NULL )
  {
  }



  template< class DiscreteFunctionType, class ProblemType >
  void modifyRHS( DiscreteFunctionType& rhs, const ProblemType& problem )
  {
    trivialProjectionOnAffineSpace( rhs, problem );
  }

  DiscreteFunctionAffinePartType& getAffinePart() const
  {
    return *df_;
  }

private:

  /**
   * @brief adds the affine part to a discrete function
   *
   * @param discreteFunction The discrete function 
   * @param problem (will be replaced)
   */
  template< class DiscreteFunctionType, class ProblemType >
  void trivialProjectionOnAffineSpace( DiscreteFunctionType& discreteFunction, const ProblemType& problem )
  {
    //TODO Lagrange-Spezialisierung ändern:
    typedef typename DiscreteFunctionType::LocalFunctionType
      LocalFunctionType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::LagrangePointSetType
      LagrangePointSetType;
    const int faceCodim = 1;
    typedef typename LagrangePointSetType::template Codim< faceCodim >::SubEntityIteratorType
      FaceDofIteratorType;
    typedef typename EntityType::Geometry
      GeometryType;
    typedef typename IntersectionIteratorType::Intersection
      IntersectionType;

    //TODO Implementation vervollstaendigen.
    ConstraintsType constraints = space_.getConstraints();

    const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();
    const GridPartType &gridPart = dfSpace.gridPart();

    bool hasDirichletBoundary = false;
    const IteratorType end = space_.end();
    for( IteratorType it = space_.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      // if entity has boundary intersections 
      if( entity.hasBoundaryIntersections() )
      {

        // get local functions of data 
        LocalFunctionType discreteFunctionLocal = discreteFunction.localFunction( entity );

        const GeometryType& geo = entity.geometry(); 

        const LagrangePointSetType &lagrangePointSet = dfSpace.lagrangePointSet( entity );

        IntersectionIteratorType it = gridPart.ibegin( entity );
        const IntersectionIteratorType endit = gridPart.iend( entity );
        for( ; it != endit; ++it )
        {
          const IntersectionType& intersection = *it;

          // if intersection is with boundary, adjust data
          if( intersection.boundary() )
          {
            hasDirichletBoundary = true;

            // get face number of boundary intersection
            const int face = intersection.indexInInside();
            // get dof iterators 
            FaceDofIteratorType faceIt
              = lagrangePointSet.template beginSubEntity< faceCodim >( face );
            const FaceDofIteratorType faceEndIt
              = lagrangePointSet.template endSubEntity< faceCodim >( face );
            for( ; faceIt != faceEndIt; ++faceIt )
            {
              // get local dof number 
              const int localDof = *faceIt;

              // get global coordinate of point on boundary
              const DomainType global = geo.global( lagrangePointSet.point( localDof ) );

              // check whether Dirichlet boundary or not
              // (remark: all boundary nodes are supposed to be Dirichlet boundary faces at the moment)
              //if( ! constraints.problem.dirichletBoundary(intersection.boundaryId(), global ) )
              // {
              //   continue;
              // }

              //! TODO: phi ersetzen durch df_ und problem entfernen

              // evaluate boundary data
              RangeType phi;
              problem.g( global, phi );
              //!const RangeFieldType xsqr = global*global;
              //!phi = exp( -10.0 * xsqr );

              // adjust right hand side and solution data
              discreteFunctionLocal[localDof] = phi[0];

            }
          }
        }
      }
    }
  }

  LinearSubspaceType& space_;
  DiscreteFunctionAffinePartType* df_;
};


/**
 * @class AffineSubspaceDiscreteFunction
 *
 * @brief A class representing a discrete function from a linear subspace.
 *
 * @tparam AffineSubspaceImp An affine subspace.
 * @tparam LinearSubspaceDiscreteFunctionImp A discrete function for a linear subspace.
 *
 * @todo improve comments
 *
 * @attention class is not finished yet
 */
template< class AffineSubspaceImp, class LinearSubspaceDiscreteFunctionImp >
class AffineSubspaceDiscreteFunction
  : public LinearSubspaceDiscreteFunctionImp
{
public:
  typedef LinearSubspaceDiscreteFunctionImp
    LinearSubspaceDiscreteFunctionType;
  typedef typename LinearSubspaceDiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef AffineSubspaceImp
    AffineSubspaceType;
 
  //TODO: assert: LinearSubspaceDiscreteFunctionType::DiscreteFunctionSpaceType and
  //              AffineSubspaceType::DiscreteFunctionSpaceType fits together
  
  

  typedef typename LinearSubspaceDiscreteFunctionType::DomainFieldType
    DomainFieldType;
  typedef typename LinearSubspaceDiscreteFunctionType::RangeFieldType
    RangeFieldType;
  typedef typename LinearSubspaceDiscreteFunctionType::DomainType
    DomainType;
  typedef typename LinearSubspaceDiscreteFunctionType::RangeType
    RangeType;

  typedef typename LinearSubspaceDiscreteFunctionType::JacobianRangeType
    JacobianRangeType;
  typedef typename LinearSubspaceDiscreteFunctionType::HessianRangeType
    HessianRangeType;

  typedef typename LinearSubspaceDiscreteFunctionType::GridType
    GridType;
  typedef typename LinearSubspaceDiscreteFunctionType::LocalFunctionType
    LocalFunctionType;
  typedef typename LinearSubspaceDiscreteFunctionType::DofIteratorType
    DofIteratorType;
  typedef typename LinearSubspaceDiscreteFunctionType::ConstDofIteratorType
    ConstDofIteratorType;
  typedef typename LinearSubspaceDiscreteFunctionType::MappingType
    MappingType;
  typedef typename LinearSubspaceDiscreteFunctionType::FunctionSpaceType
    FunctionSpaceType;
  typedef typename LinearSubspaceDiscreteFunctionType::FunctionType
    FunctionType;

private:
  typedef LinearSubspaceDiscreteFunctionImp
    BaseType;

public:

  /**
   * @brief constructor
   *
   * The constructor automatically adds the affine part described by
   * the AffineSubspace
   *
   * @param affDiscFuncSpace discrete function space (here: affine subspace)
   * @param linDiscFunc discrete function (here: discrete function for linear subspace)
   */
  AffineSubspaceDiscreteFunction( AffineSubspaceType& affDiscFuncSpace, 
                                  LinearSubspaceDiscreteFunctionType& linDiscFunc)
    : affDiscFuncSpace_( affDiscFuncSpace ), linDiscFunc_( linDiscFunc )
  {
    //TODO problem muss noch raus
    //discFuncSpace_.modifyRHS( discFunc_, problem );
  }

  /**
   * @brief Returns the discrete function.
   */
  LinearSubspaceDiscreteFunctionType& discreteFunction()
  {
    return linDiscFunc_;
  }

  /**
   * @brief Returns the discrete function space (here: linear subspace).
   */

  AffineSubspaceType& affineSubspace()
  {
    return affDiscFuncSpace_;
  }

private:
  AffineSubspaceType& affDiscFuncSpace_;
  LinearSubspaceDiscreteFunctionType& linDiscFunc_;

};

// A small conceptual view of how the function spaces/discrete functions are created:
// Normally a DiscreteFunction is created depending on a DiscreteFunctionSpace.
// The subspaces are created from existing function spaces, discrete function analog.
//
// DiscreteFunction                 <----     DiscreteFunctionSpace
//
//          |                                         | (Constraints)
//          v                                         v
//
// LinearSubspaceDiscreteFunction   <----     LinearSubspace
//
//          |                                         | (DiscreteFunctionAffinePart)
//          v                                         v
//
// AffineSubspaceDiscreteFunction   <----     AffineSubspace
//
//
// Every class has got a public typedef DiscreteFunctionSpaceType 
// ("pointing" to the basic DiscreteFunctionSpace).
//
// I hope it's clear by now...
//



template< class GridImp, int polynomialOrder,
          class GridPartImp = Dune::AdaptiveLeafGridPart< GridImp, Dune::InteriorBorder_Partition > >
class Algorithm                                                       /*@LST0S@*/
{
public:                                                               /*@LST0E@*/
  typedef GridImp
    HGridType;

  /** choose the grid partition (and hence the index set) to use
   *
   *  \note Not all index sets are continuous. The LeafIndexSet for AlbertaGrid,
   *        for example, is not. If you want to use OEM solvers, the index set
   *        must be continuous. In such a case use AdaptiveLeafGridPart.
   */
  //---- GridParts -----------------------------------------------------------
  typedef GridPartImp
    GridPartType;

  //---- FunctionSpace -------------------------------------------------------
  //! define the function space, \f[ \R^n \rightarrow \R \f]
  // see dune/common/functionspace.hh
  typedef Dune::FunctionSpace< double, double, HGridType::dimensionworld, 1 >
    FunctionSpaceType;

  //! to be revised 
  typedef Dune::ProblemInterface< FunctionSpaceType >
    ProblemType;

  // The data functions (as defined in problem.hh)
  //---- Right Hand Side, Exact Solution, and Stiffness tensor ---------------
  typedef typename ProblemType::ExactSolutionType
    ExactSolutionType;

  //---- Adapter for exact solution ------------------------------------------
  typedef Dune::DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
    GridExactSolutionType;

  //---- DiscreteFunctionSpace -----------------------------------------------
  //! define the discrete function space our unkown belongs to
  typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, 
                                               GridPartType, 
                                               polynomialOrder, 
                                               Dune::CachingStorage >
    DiscreteSpaceType;

  //---- Functional ---------------------------------------------------
  //! define the functionals you want to use
  //typedef Dune::LinearCodimZeroFunctional< FunctionSpaceType, DiscreteSpaceType >
  //  LinearFunctionalType;


  //----- Constraints ----------------------------------------------------------
  //! define the Constraints you want to use, i.e. Dirichlet constraints
  //typedef VectorConstraints< LinearFunctionalType >
  //  ConstraintsType;
  typedef Dune::DirichletConstraints< DiscreteSpaceType >
    DirichletConstraintsType;

  //---- DiscreteFunction ----------------------------------------------------
  //---- good choice for adaptive simulations using OEM solver ---------------
  //! define the type of discrete function we are using
  typedef Dune::AdaptiveDiscreteFunction< DiscreteSpaceType >
    DiscreteFunctionType;
  //---- other possible choices, use BlockVectorDiscreteFunction for ISTL ----
  //typedef Dune::BlockVectorDiscreteFunction< DiscreteSpaceType >
  //  DiscreteFunctionType;
  //typedef Dune::ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteSpaceType, DynamicVector< double > > >
  //  DiscreteFunctionType;

  //----- Subspaces ------------------------------------------------------------------- 
  typedef LinearSubspace< DiscreteSpaceType, DirichletConstraintsType >
    LinearSubspaceType;
  typedef AffineSubspace< LinearSubspaceType, DiscreteFunctionType >
    AffineSubspaceType;


  //---- MatrixObjects -------------------------------------------------------
  //---- good choice for build in solvers ------------------------------------
  //! define the type of the system matrix object
  typedef Dune::SparseRowMatrixTraits< DiscreteSpaceType, DiscreteSpaceType >
    MatrixObjectTraits;

  //---- other choices, ISTLMatrixTraits for BCRSMatrix from DUNE-ISTL -------
  //typedef Dune::ISTLMatrixTraits < DiscreteSpaceType, DiscreteSpaceType >
  //  MatrixObjectTraits;
  //typedef Dune::OnTheFlyMatrixTraits < DiscreteSpaceType, DiscreteSpaceType >
  //  MatrixObjectTraits;

  //! define the discrete laplace operator, see ./laplaceoperator.hh 
  typedef Dune::DifferentialOperator< LinearSubspaceType, AffineSubspaceType, MatrixObjectTraits >
    DifferentialOperatorType;


  //---- InverseOperator ----------------------------------------------------
  //---- good choice for build in CG solver ---------------------------------
  //! define the inverse operator we are using to solve the system
  typedef Dune::CGInverseOp< DiscreteFunctionType, DifferentialOperatorType >
    InverseOperatorType;
    //---- other choices ------------------------------------------------------
    //typedef Dune::ISTLBICGSTABOp< DiscreteFunctionType, DifferentialOperatorType >
    //typedef Dune::ISTLCGOp< DiscreteFunctionType, DifferentialOperatorType >
    //typedef Dune::OEMCGOp< DiscreteFunctionType,DifferentialOperatorType >
    //typedef Dune::OEMBICGSTABOp< DiscreteFunctionType,DifferentialOperatorType >
    //typedef Dune::OEMBICGSQOp< DiscreteFunctionType,DifferentialOperatorType >
    //typedef Dune::OEMGMRESOp< DiscreteFunctionType,DifferentialOperatorType >
    //  InverseOperatorType;

public:

  //! constructor                                                  /*@LST0S@*/
  Algorithm( HGridType& grid, const ProblemType& problem )
    : grid_( grid ),
      gridPart_( grid_ ),
      space_( gridPart_ ),
      problem_( problem ),
      constraints_( space_ ),
      linSubspace_( space_, constraints_ ),
      affSubspace_( linSubspace_ )
  {
    // add entries to eoc calculation
    std::vector< std::string > femEocHeaders;
    // we want to calculate L2 and H1 error (and EOC)
    femEocHeaders.push_back( "$L^2$-error" );
    femEocHeaders.push_back( "$H^1$-error" );                        /*@LST0E@*/

    // get eoc id
    eocId_ = Dune::FemEoc::addEntry( femEocHeaders );

    // distribute grid (only for parallel runs, in serial runs nothing is done here)
    grid.loadBalance();
  }

  //! setup and solve the linear system
  void operator()( DiscreteFunctionType& solution )
  {                                                                /*@LST0E@*/
    // functional describing the right hand side
    typedef Dune::IntegralFunctional< ProblemType, DiscreteFunctionType >
      FunctionalType;

    // in verbose mode some info 
    if( Dune::Parameter::verbose() ) 
    {
      std::cout << std::endl << "Solving for " << space_.size()
        << " unkowns and polynomial order "
        << DiscreteSpaceType::polynomialOrder << "."
        << std::endl << std::endl;
    }

    // initialize solution with zero                                 /*@LST0S@*/
    solution.clear();

    // create laplace assembler (is assembled on call of systemMatrix by solver)
    DifferentialOperatorType laplace( linSubspace_ , affSubspace_, problem_ );       /*@\label{poi:laplaceDef}@*/

    FunctionalType rhsFunctional( problem_ , space_ );

    // solve the linear system  @\ref{
    solve( laplace, rhsFunctional, solution );  /*@\label{poi:solve}@*/
  }

  //! finalize computation by calculating errors and EOCs
  void finalize( DiscreteFunctionType& solution )
  {
    // create exact solution
    ExactSolutionType uexact( problem_ );

    // create grid function adapter
    GridExactSolutionType ugrid( "exact solution", uexact, gridPart_,
                                 DiscreteSpaceType::polynomialOrder + 1 );

    // create L2 - Norm
    Dune::L2Norm< GridPartType > l2norm( gridPart_ );
    // calculate L2 - Norm
    const double l2error = l2norm.distance( ugrid, solution );

    // create H1 - Norm
    Dune::H1Norm< GridPartType > h1norm( gridPart_ );
    // calculate H1 - Norm
    const double h1error = h1norm.distance( ugrid, solution );

    // store values
    std::vector< double > errors;
    errors.push_back( l2error );
    errors.push_back( h1error );

    // submit error to the FEM EOC calculator
    Dune::FemEoc::setErrors( eocId_, errors );
  }

  //! return reference to discrete space 
  DiscreteSpaceType& space()
  {
    return space_;
  }                    /*@LST0E@*/

private:                                                            /*@LST0S@*/
  //! solve the resulting linear system
  template< class FunctionalType >
  void solve( DifferentialOperatorType& laplace,
              const FunctionalType& functional,
              DiscreteFunctionType& solution )
  {                                                                /*@LST0E@*/
    // create storage for rhs (this could be any vector type)
    DiscreteFunctionType rhs( "rhs", space_ );
    rhs.clear();

    // assemble right hand side functional vector
    functional.algebraic( rhs );
    affSubspace_.modifyRHS( rhs, problem_ );

    //TODO das hier muss geändert werden (solution aus AffineDiscreteSubspace 
    //hat damit automatisch die richtigen Randwerte)
    affSubspace_.modifyRHS( solution, problem_ );

    // create timer (also stops time)
    Dune::Timer timer;

    // solve the linear system (with CG)
    const double reduction = 1e-6;
    double solverEps = 1e-8;
    solverEps = Dune::Parameter::getValue( "femhowto.solvereps", solverEps );

    // get verbosity information from Parameter class
    const bool verbose = Dune::Parameter::verbose();

    // create inverse operator (note reduction is not used here)    /*@LST0S@*/
    InverseOperatorType cg( laplace, reduction, solverEps );/*@\label{poi:cgInit}@*/

    // solve the system
    cg( rhs, solution );                                    /*@\label{poi:cgExec}@*//*@LST0E@*/

    // get elapsed time
    const double solveTime = timer.elapsed();

    // output in only in verbose mode (parameter fem.verbose)
    if( verbose )
      std::cout << "Time needed by solver: " << solveTime << "s" << std::endl;
  }                                                                /*@LST0S@*/

protected:
  HGridType& grid_;            // reference to grid, i.e. the hierarchical grid 
  GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid 
  DiscreteSpaceType space_;   // the discrete function space
  const ProblemType& problem_; // the problem data
  const DirichletConstraintsType constraints_;
  LinearSubspaceType linSubspace_;
  AffineSubspaceType affSubspace_;
  int eocId_;                  // id for FemEOC
};                                                                /*@LST0E@*/

#endif // POISSONSERIAL_HH
