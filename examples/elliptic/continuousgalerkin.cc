#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/function/expression.hh>

#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/linear.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed/discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed/discretizations/la/factory/eigen.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed/discretizations/assembler/system/constrained.hh>
#include <dune/detailed/discretizations/la/backend/solver/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>

const std::string id = "elliptic.continuousgalerkin";

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
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << ".grid" << std::endl;
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "numElements = 4" << std::endl;
    file << "[diffusion]" << std::endl;
    file << "order = 0"  << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[force]" << std::endl;
    file << "order = 0"  << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[dirichlet]" << std::endl;
    file << "order = 0"  << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [0.0; 0.0; 0.0]" << std::endl;
    file << "[solver]" << std::endl;
    file << "maxIter = 5000"  << std::endl;
    file << "precision = 1e-12"  << std::endl;
    file << "[visualization]" << std::endl;
    file << "filename = " << id << ".solution" << std::endl;
    file << "name = solution" << std::endl;
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
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".param";
    ensureParamFile(paramFilename);
    Dune::Stuff::Common::ExtendedParameterTree paramTree(argc, argv, paramFilename);
    if (!paramTree.hasSub(id))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing sub " << id << " in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));

    // timer
    Dune::Timer timer;

    // grid
    std::cout << "setting up grid:" << std::endl;
    typedef Dune::Stuff::Grid::Provider::Cube<> GridProviderType;
    const GridProviderType gridProvider = GridProviderType::createFromParamTree(paramTree);
    typedef GridProviderType::GridType GridType;
    const Dune::shared_ptr< const GridType > grid = gridProvider.grid();
    typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
    const GridPartType gridPart(*grid);
    std::cout << "  took " << timer.elapsed() << " sec, has " << grid->size(0) << " entities" << std::endl;
    std::cout << "visualizing grid... " << std::flush;
    if (!paramTree.sub(id).hasKey("filename"))
      DUNE_THROW(Dune::RangeError,
                 "\nError: missing key 'filename' in the following Dune::parameterTree:\n" << paramTree.sub(id).reportString("  "););
    timer.reset();
    gridProvider.visualize(paramTree.sub(GridProviderType::id()).get("filename", id + ".grid"));
    std::cout << " done (took " << timer.elapsed() << " sek)" << std::endl;

    // spaces
    std::cout << "initializing spaces... " << std::flush;
    timer.reset();
    const int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const int DUNE_UNUSED(dimRange) = 1;
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType, GridPartType, polOrder > TestSpaceType;
    const TestSpaceType testSpace(gridPart);
    typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Sub::Linear::Dirichlet< TestSpaceType > AnsatzSpaceType;
    const AnsatzSpaceType ansatzSpace(testSpace);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // data
    std::cout << "initializing data functions... " << std::flush;
    timer.reset();
    typedef Dune::Stuff::Function::Expression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ExpressionFunctionType;
    if (!paramTree.hasSub("diffusion"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing sub 'diffusion' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
    if (!paramTree.sub("diffusion").hasKey("order"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'order' in the following Dune::ParameterTree:\n" << paramTree.sub("order").reportString("  "));
    const Dune::shared_ptr< const ExpressionFunctionType >
        diffusion(new ExpressionFunctionType(ExpressionFunctionType::createFromParamTree(paramTree.sub("diffusion"))));
    if (!paramTree.hasSub("force"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing sub 'force' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
    if (!paramTree.sub("force").hasKey("order"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'order' in the following Dune::ParameterTree:\n" << paramTree.sub("force").reportString("  "));
    const Dune::shared_ptr< const ExpressionFunctionType >
        force(new ExpressionFunctionType(ExpressionFunctionType::createFromParamTree(paramTree.sub("force"))));
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // left hand side (operator)
    std::cout << "initializing operator and functional... " << std::flush;
    timer.reset();
    typedef Dune::Detailed::Discretizations
        ::Evaluation
        ::Local
        ::Binary
        ::Elliptic< FunctionSpaceType, ExpressionFunctionType >
      EllipticEvaluationType;
    const EllipticEvaluationType ellipticEvaluation(diffusion, paramTree.sub("diffusion").get("order", 0));
    typedef Dune::Detailed::Discretizations
        ::DiscreteOperator
        ::Local
        ::Codim0
        ::Integral< EllipticEvaluationType >
      EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(ellipticEvaluation);
    // right hand side (functional)
    typedef Dune::Detailed::Discretizations
        ::Evaluation
        ::Local
        ::Unary
        ::Scale< FunctionSpaceType, ExpressionFunctionType >
      ProductEvaluationType;
    const ProductEvaluationType productEvaluation(force, paramTree.sub("force").get("order", 0));
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunctional
        ::Local
        ::Codim0
        ::Integral< ProductEvaluationType >
      L2FunctionalType;
    const L2FunctionalType l2Functional(productEvaluation);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // system matrix and right hand side
    std::cout << "initializing matrix and vector containers..." << std::flush;
    timer.reset();
    typedef AnsatzSpaceType::PatternType PatternType;
    const Dune::shared_ptr< const PatternType > pattern = ansatzSpace.computePattern(testSpace);
    typedef Dune::Detailed::Discretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;
    typedef ContainerFactory::SparseMatrixType MatrixBackendType;
    MatrixBackendType systemMatrix = ContainerFactory::createSparseMatrix(ansatzSpace.map().size(), testSpace.map().size(), *pattern);
    typedef ContainerFactory::DenseVectorType VectorBackendType;
    VectorBackendType rhs = ContainerFactory::createDenseVector(testSpace.map().size());
    VectorBackendType ret = ContainerFactory::createDenseVector(testSpace);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // assembler
    std::cout << "assembing system... " << std::flush;
    timer.reset();
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
    const LocalMatrixAssemblerType localmatrixAssembler(ellipticOperator);
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2FunctionalType > LocalVectorAssemblerType;
    const LocalVectorAssemblerType localVectorAssembler(l2Functional);
    typedef Dune::Detailed::Discretizations::Assembler::System::Constrained< AnsatzSpaceType, TestSpaceType > SystemAssemblerType;
    const SystemAssemblerType systemAssembler(ansatzSpace, testSpace);
    systemAssembler.assembleSystem(localmatrixAssembler, systemMatrix, localVectorAssembler, rhs);
    systemAssembler.applyConstraints(systemMatrix, rhs);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solve system
    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabIlut Solver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::BicgstabDiagonal Solver;
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalUpper Solver; // seems to produce strange results
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::CgDiagonalLower Solver; // seems to produce strange results
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyUpper Solver; // seems to produce strange results
//    typedef Dune::Detailed::Discretizations::LA::Solver::Eigen::SimplicialcholeskyLower Solver; // seems to produce strange results
    std::cout << "solving linear system (of size " << systemMatrix.rows() << "x" << systemMatrix.cols() << ")" << std::endl
              << "  using " << Solver::id << "... " << std::flush;
    if (!paramTree.hasSub("solver"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing sub 'solver' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
    if (!paramTree.sub("solver").hasKey("maxIter"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'maxIter' in the following Dune::ParameterTree:\n" << paramTree.sub("solver").reportString("  "));
    if (!paramTree.sub("solver").hasKey("precision"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'precision' in the following Dune::ParameterTree:\n" << paramTree.sub("solver").reportString("  "));
    timer.reset();
    Solver::apply(
      systemMatrix,
      ret,
      rhs,
      paramTree.sub("solver").get("maxIter", 5000),
      paramTree.sub("solver").get("precision", 1e-12));
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // postprocess
    if (!paramTree.hasSub("visualization"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing sub 'visualization' in the following Dune::ParameterTree:\n" << paramTree.reportString("  "));
    if (!paramTree.sub("visualization").hasKey("name"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'name' in the following Dune::ParameterTree:\n" << paramTree.sub("visualization").reportString("  "));
    if (!paramTree.sub("visualization").hasKey("filename"))
      DUNE_THROW(Dune::RangeError,
                 "\Error: missing key 'filename' in the following Dune::ParameterTree:\n" << paramTree.sub("visualization").reportString("  "));
    const std::string solutionName = paramTree.sub("visualization").get("name", "solution");
    const std::string solutionFilename = paramTree.sub("visualization").get("filename", id + ".solution");
    std::cout << "writing '" << solutionName << "' to '" << solutionFilename;
    if (dimDomain == 1)
      std::cout << ".vtp";
    else
      std::cout << ".vtu";
    std::cout << "'... " << std::flush;
    typedef Dune::Detailed::Discretizations::DiscreteFunction::Default< AnsatzSpaceType, VectorBackendType > DiscreteFunctionType;
    Dune::shared_ptr< DiscreteFunctionType > solution(new DiscreteFunctionType(ansatzSpace, ret, solutionName));
    typedef Dune::VTKWriter< AnsatzSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(ansatzSpace.gridView());
    vtkWriter.addVertexData(solution);
    vtkWriter.write(solutionFilename);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // done
    return 0;
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
