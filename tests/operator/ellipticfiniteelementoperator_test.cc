
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

// system includes
#include <iostream>

// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>

// dune common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune fem includes
//#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

// reenable warnings about problems in dune headers
#include <dune/fem-tools/header/enablewarnings.hh>

// dune fem-functionals includes
#include <dune/fem/operator/ellipticfiniteelement.hh>

// dune fem-tools includes
#include <dune/fem-tools/function/functiontools.hh>

/**
  * \brief Analytical function which induces the functional.
  **/
template < class FunctionSpaceImp >
class AnalyticalFunction : public Dune::Fem::Function< FunctionSpaceImp , AnalyticalFunction < FunctionSpaceImp > >
{
  public:
    typedef FunctionSpaceImp
      FunctionSpaceType;
    typedef AnalyticalFunction< FunctionSpaceType >
      ThisType;
    typedef Dune::Function< FunctionSpaceType, ThisType >
      BaseType;
    typedef typename FunctionSpaceType::DomainType
      DomainType;
    typedef typename FunctionSpaceType::RangeType
      RangeType;
    typedef typename FunctionSpaceType::RangeFieldType
      RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType
      JacobianRangeType;

    AnalyticalFunction(){}

    ~AnalyticalFunction(){}

    inline void evaluate( const DomainType& arg, RangeType& ret ) const
    {
      ret = 1.0;
    }

    inline void jacobian( const DomainType& arg, JacobianRangeType& ret ) const
    {
      ret = JacobianRangeType( 1.0 );
    }

};

// disable warnings about problems in dune headers (sourced by dgfparser)
#include <dune/fem-tools/header/disablewarnings.hh>

// main
int main(int argc, char** argv)
{
  try{

    // print welcome
    std::cout << "Elliptic finite element operator test ";

    // mpi
    Dune::MPIManager::initialize ( argc, argv );

    // dimension and grid
    const int dim = 2;

    typedef Dune::YaspGrid< dim >
      GridType;

    Dune::GridPtr< GridType > gridPtr( "unitcube_2d.dgf" );

    typedef Dune::AdaptiveLeafGridPart< GridType >
      GridPartType;

    GridPartType gridPart( *gridPtr );

    // analytical function space and function
    typedef Dune::FunctionSpace< double, double, dim, 1 >
      AnalyticalFunctionSpaceType;
    typedef AnalyticalFunction< AnalyticalFunctionSpaceType >
      AnalyticalFunctionType;

    const AnalyticalFunctionType analyticalFunction;

    // discrete function space and function
    const int polOrder = 1;

    typedef Dune::LagrangeDiscreteFunctionSpace< AnalyticalFunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;

    const DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

    typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
      DiscreteFunctionType;

    DiscreteFunctionType discreteFunction( "discrete_function" , discreteFunctionSpace );
    Dune::FemTools::setDiscreteFunctionToScalarValue( discreteFunction, 1.0 );

    // operator
    typedef Dune::Functionals::Operator::EllipticFiniteElement< DiscreteFunctionSpaceType, AnalyticalFunctionType >
      EllipticFiniteElementOperatorType;

    EllipticFiniteElementOperatorType ellipticFiniteElementOperator( discreteFunctionSpace, analyticalFunction);

    // entity
    EllipticFiniteElementOperatorType::EntityIteratorType entityIterator = discreteFunctionSpace.begin();
    EllipticFiniteElementOperatorType::EntityType& entity = *entityIterator;

    // test applyLocal
    ellipticFiniteElementOperator.applyLocal( entity );

//    // test as scalar product
//    const double volume = ellipticFiniteElementOperator( discreteFunction, discreteFunction );

//    // functions are chosen to equal 1 when multiplied, thus the application of the operator should yield the volume
//    // of the area, which in turn should be 1 in case of the two-dimensional unitcube
//    if ( volume == 1.0 )
//      std::cout << "passed!" << std::endl;
//    else
//      std::cout << "failed (result should equal 1, is " << volume << ")!" << std::endl;

    // we don't make no errors^^
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
} // end main
