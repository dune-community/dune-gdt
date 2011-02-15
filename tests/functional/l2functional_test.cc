#ifdef HAVE_CONFIG_H
#include "config.h"
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
#include <dune/fem/functional/l2functional.hh>

// dune fem-tools includes
#include "../../tools/function/functiontools.hh"

/**
  * \brief Analytical function which induces the functional.
  **/
template <class FunctionSpaceImp>
class AnalyticalFunction : public Dune::Function<FunctionSpaceImp, AnalyticalFunction<FunctionSpaceImp>>
{
public:
  typedef AnalyticalFunction<FunctionSpaceImp> ThisType;
  typedef Dune::Function<FunctionSpaceImp, ThisType> BaseType;
  //    typedef typename BaseType::DomainType
  //      DomainType;
  //    typedef typename BaseType::RangeType
  //      RangeType;
  typedef typename FunctionSpaceImp::RangeFieldType RangeFieldType;

  AnalyticalFunction()
  {
  }

  ~AnalyticalFunction()
  {
  }

  template <class DomainType, class RangeType>
  inline void evaluate(const DomainType& arg, RangeType& ret) const
  {
    ret = 1.0;
  }
};

// main
int main(int argc, char** argv)
{
  try {

    // print welcome
    std::cout << "This is l2functional_test." << std::endl;

    // mpi
    Dune::MPIManager::initialize(argc, argv);

    //    // command line arguments
    //    Dune::Parameter::append( argc, argv );

    // dimension and grid
    const int dim = 2;

    typedef Dune::YaspGrid<dim> GridType;

    Dune::GridPtr<GridType> gridPtr("unitcube_2d.dgf");

    typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;

    GridPartType gridPart(*gridPtr);

    // analytical function space and function
    typedef Dune::FunctionSpace<double, double, dim, 1> AnalyticalFunctionSpaceType;
    typedef AnalyticalFunction<AnalyticalFunctionSpaceType> AnalyticalFunctionType;

    const AnalyticalFunctionType analyticalFunction;

    // discrete function space and function
    const int polOrder = 1;

    typedef Dune::LagrangeDiscreteFunctionSpace<AnalyticalFunctionSpaceType, GridPartType, polOrder>
        DiscreteFunctionSpaceType;

    const DiscreteFunctionSpaceType discreteFunctionSpace(gridPart);

    typedef Dune::AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;

    DiscreteFunctionType discreteFunction("discrete_function", discreteFunctionSpace);
    Dune::FemTools::setDiscreteFunctionToScalarValue(discreteFunction, 1.0);

    // test functional
    std::cout << "Testing L2 functional: ";

    typedef Dune::Functionals::L2Functional<AnalyticalFunctionType> L2FunctionalType;
    const L2FunctionalType l2Functional(analyticalFunction);

    const double volume = l2Functional(discreteFunction);

    if (volume == 1.0)
      std::cout << " test passed!" << std::endl;
    else
      std::cout << " test failed (result should equal 1, is " << volume << ")!" << std::endl;

    // we don't make no errors^^
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
