
#include "config.h"

// system
#include <iostream>
#include <sstream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/helper-tools/header/disablewarnings.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/helper-tools/header/enablewarnings.hh>

// dune-helper-tools
#include <dune/helper-tools/common/parametertree.hh>
#include <dune/helper-tools/grid/provider/cube.hh>
#include <dune/helper-tools/function/expression.hh>

// dune-detailed-discretizations
#include <dune/detailed-discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed-discretizations/discretefunctionspace/subspace/linear.hh>
#include <dune/detailed-discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim0/integral.hh>

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[helper-tools.grid.provider.cube]" << std::endl;
    file << "level = 2" << std::endl;
    file << "[data.a]" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression.0 = 1.0" << std::endl;
    file << "expression.1 = 1.0" << std::endl;
    file << "expression.2 = 1.0" << std::endl;
    file << "order = 0" << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);
    // parameter
    const std::string filename = "continuous_galerkin.param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::HelperTools::Common::ParameterTree::init(argc, argv, filename);
    // grid
    typedef Dune::HelperTools::Grid::Provider::UnitCube<Dune::GridSelector::GridType> GridProviderType;
    GridProviderType gridProvider(paramTree);
    typedef GridProviderType::GridType GridType;
    GridType& grid = gridProvider.grid();
    typedef Dune::LeafGridPart<GridType> GridPartType;
    GridPartType gridPart(grid);
    // function spaces
    const int dimDomain = GridProviderType::dim;
    const int dimRange  = 1;
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    typedef Dune::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, dimRange> FunctionSpaceType;
    typedef Dune::DetailedDiscretizations::DiscreteFunctionSpace::Continuous::Lagrange<FunctionSpaceType,
                                                                                       GridPartType,
                                                                                       polOrder> DiscreteH1Type;
    const DiscreteH1Type discreteH1(gridPart);
    typedef Dune::DetailedDiscretizations::DiscreteFunctionSpace::Subspace::Linear::Dirichlet<DiscreteH1Type>
        AnsatzSpaceType;
    const AnsatzSpaceType ansatzSpace(discreteH1);
    typedef AnsatzSpaceType TestSpaceType;
    const TestSpaceType testSpace(discreteH1);
    // left hand side operator hand side
    typedef Dune::HelperTools::Function::Expression<DomainFieldType, dimDomain, RangeFieldType, dimRange>
        DataFunctionType;
    if (!paramTree.hasSub("data.a")) {
      std::stringstream msg;
      msg << "Error in " << filename << ": key 'data.a' not found in the following Dune::Parametertree" << std::endl;
      paramTree.report(msg);
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    const DataFunctionType a(paramTree.sub("data.a"));
    typedef Dune::DetailedDiscretizations::Evaluation::Local::Binary::Elliptic<FunctionSpaceType>
        EllipticEvaluationType;
    EllipticEvaluationType ellipticEvaluation(a);
    typedef Dune::DetailedDiscretizations::DiscreteOperator::Local::Codim0::Integral<EllipticEvaluationType>
        EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(ellipticEvaluation);


    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
