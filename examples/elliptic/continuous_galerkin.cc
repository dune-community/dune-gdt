
#include "config.h"

// system
#include <iostream>
#include <sstream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

// dune-grid
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// dune-fem
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/gridpartview.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/function/expression.hh>

// dune-detailed-discretizations
#include <dune/detailed-discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed-discretizations/discretefunctionspace/subspace/linear.hh>
#include <dune/detailed-discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed-discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed-discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed-discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed-discretizations/la/factory/eigen.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed-discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed-discretizations/assembler/system/constrained.hh>
#include <dune/detailed-discretizations/la/backend/solver/eigen.hh>
#include <dune/detailed-discretizations/discretefunction/default.hh>

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
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "level = 4" << std::endl;
    file << "[diffusion]" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression.0 = 1.0"  << std::endl;
    file << "expression.1 = 1.0"  << std::endl;
    file << "expression.2 = 1.0"  << std::endl;
    file << "order = 0"  << std::endl;
    file << "[force]" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression.0 = 1.0"  << std::endl;
    file << "expression.1 = 1.0"  << std::endl;
    file << "expression.2 = 1.0"  << std::endl;
    file << "order = 0"  << std::endl;
    file << "[solver]" << std::endl;
    file << "maxIter = 5000"  << std::endl;
    file << "precision = 1e-12"  << std::endl;
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
    const std::string id = "continuous_galerkin";
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::Stuff::Common::Parameter::Tree::init(argc, argv, filename);

    // timer
    Dune::Timer timer;

    // grid
    std::cout << "setting up grid:" << std::endl;
    typedef Dune::Stuff::Grid::Provider::UnitCube< Dune::GridSelector::GridType > GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, id);
    GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::GridType GridType;
    GridType& grid = gridProvider.grid();
    typedef Dune::LeafGridPart< GridType > GridPartType;
    GridPartType gridPart(grid);
    typedef Dune::GridPartView< GridPartType > GridViewType;
    GridViewType gridView(gridPart);
    std::cout << "took " << timer.elapsed() << " sec, has " << gridView.size(0) << " entities" << std::endl;

    // function spaces
    std::cout << "setting up function spaces... " << std::flush;
    timer.reset();
    const int dimDomain = GridProviderType::dim;
    const int dimRange = 1;
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    typedef Dune::DetailedDiscretizations::DiscreteFunctionSpace::Continuous::Lagrange< FunctionSpaceType, GridPartType, polOrder > DiscreteH1Type;
    const DiscreteH1Type discreteH1(gridPart);
    typedef Dune::DetailedDiscretizations::DiscreteFunctionSpace::Subspace::Linear::Dirichlet< DiscreteH1Type > AnsatzSpaceType;
    const AnsatzSpaceType ansatzSpace(discreteH1);
    typedef AnsatzSpaceType TestSpaceType;
    const TestSpaceType testSpace(discreteH1);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // left hand side (operator)
    std::cout << "setting up operator and functional... " << std::flush;
    timer.reset();
    typedef Dune::DetailedDiscretizations::Evaluation::Local::Binary::Elliptic< FunctionSpaceType > EllipticEvaluationType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "diffusion", id);
    const EllipticEvaluationType ellipticEvaluation(paramTree.sub("diffusion"));
    typedef Dune::DetailedDiscretizations::DiscreteOperator::Local::Codim0::Integral< EllipticEvaluationType > EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(ellipticEvaluation);

    // right hand side (functional)
    typedef Dune::DetailedDiscretizations::Evaluation::Local::Unary::Scale< FunctionSpaceType > ProductEvaluationType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "force", id);
    const ProductEvaluationType productEvaluation(paramTree.sub("force"));
    typedef Dune::DetailedDiscretizations::DiscreteFunctional::Local::Codim0::Integral< ProductEvaluationType > L2FunctionalType;
    const L2FunctionalType l2Functional(productEvaluation);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // system matrix and right hand side
    std::cout << "setting up matrix and vector container... " << std::flush;
    timer.reset();
    typedef Dune::DetailedDiscretizations::LA::Factory::Eigen< RangeFieldType > ContainerFactory;
    typedef ContainerFactory::SparseMatrixType MatrixType;
    MatrixType systemMatrix = ContainerFactory::createSparseMatrix(ansatzSpace, testSpace);
    typedef ContainerFactory::DenseVectorType VectorType;
    VectorType rhs = ContainerFactory::createDenseVector(testSpace);
    VectorType solution = ContainerFactory::createDenseVector(testSpace);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // assembler
    std::cout << "setting up assembler... " << std::flush;
    timer.reset();
    typedef Dune::DetailedDiscretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType > LocalMatrixAssemblerType;
    const LocalMatrixAssemblerType localmatrixAssembler(ellipticOperator);
    typedef Dune::DetailedDiscretizations::Assembler::Local::Codim0::Vector< L2FunctionalType > LocalVectorAssemblerType;
    const LocalVectorAssemblerType localVectorAssembler(l2Functional);
    typedef Dune::DetailedDiscretizations::Assembler::System::Constrained< AnsatzSpaceType, TestSpaceType > SystemAssemblerType;
    const SystemAssemblerType systemAssembler(ansatzSpace, testSpace);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // assemble system
    std::cout << "assembling system... " << std::flush;
    timer.reset();
    systemAssembler.assembleSystem(localmatrixAssembler, systemMatrix, localVectorAssembler, rhs);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // solve system
//    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::BicgstabIlut Solver;
//    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::BicgstabDiagonal Solver;
    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::CgDiagonalUpper Solver;
//    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::CgDiagonalLower Solver;
//    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::SimplicialcholeskyUpper Solver;
//    typedef Dune::DetailedDiscretizations::LA::Solver::Eigen::SimplicialcholeskyLower Solver;
    std::cout << "solving linear system using " << Solver::id << "... " << std::flush;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, "solver", id);
    timer.reset();
    Solver::apply(
      systemMatrix,
      solution,
      rhs,
      paramTree.sub("solver").get("maxIter", 5000),
      paramTree.sub("solver").get("precision", 1e-12));
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // postprocess
    typedef Dune::DetailedDiscretizations::DiscreteFunction::Default< AnsatzSpaceType, VectorType > DiscreteFunctionType;
    Dune::shared_ptr< DiscreteFunctionType > u(new DiscreteFunctionType(ansatzSpace, solution, "solution"));
    typedef Dune::VTKWriter< AnsatzSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(ansatzSpace.gridView());
    vtkWriter.addVertexData(u);
    vtkWriter.write(id + "_solution");

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
