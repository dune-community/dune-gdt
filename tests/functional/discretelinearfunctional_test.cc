#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

// system includes
#include <iostream>

// dune common includes
#include <dune/common/exceptions.hh>

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

// dune fem-functionals includes
#include <dune/fem/functional/discretelinearfunctional.hh>

// dune fem-tools includes
#include "../../tools/function/functiontools.hh"

/**
  * \brief Analytical function which induces the functional.
  **/
template < class FunctionSpaceImp >
class AnalyticalFunction : public Dune::Function < FunctionSpaceImp , AnalyticalFunction < FunctionSpaceImp > >
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

    AnalyticalFunction(){}

    ~AnalyticalFunction(){}

    inline void evaluate( const DomainType& arg, RangeType& ret ) const
    {
      ret = 0.5;
    }

};

// main
int main(int argc, char** argv)
{
  try{

    // print welcome
    std::cout << "Discrete linear functional test:" << std::endl;

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
    Dune::FemTools::setDiscreteFunctionToScalarValue( discreteFunction, 2.0 );

    // test functional
//    typedef Dune::Functionals::DiscreteLinearFunctionalDefaultTraits< AnalyticalFunctionType >
//      DiscreteLinearFunctionalDefaultTraitsType;
//    typedef Dune::Functionals::DiscreteLinearFunctionalDefault< DiscreteLinearFunctionalDefaultTraitsType >
//      DiscreteLinearFunctionalDefaultType;
    typedef Dune::Functionals::TestLinearFunctionalTraits< AnalyticalFunctionType >
      TestLinearFunctionalTraitsType;
    typedef Dune::Functionals::TestLinearFunctional< TestLinearFunctionalTraitsType >
      TestLinearFunctionalType;

//    DiscreteLinearFunctionalDefaultType discreteLinearFunctionalDefault;

//    discreteLinearFunctionalDefault( discreteFunction );

    TestLinearFunctionalType testLinearFunctional( analyticalFunction );

    testLinearFunctional( discreteFunction );

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
}
