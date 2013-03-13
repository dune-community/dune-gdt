#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#define HAVE_DUNE_DETAILED_DISCRETIZATIONS 1

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
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/linear.hh>
#include <dune/detailed/discretizations/la/container/factory/eigen.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>
#include <dune/detailed/discretizations/discretefunctionspace/sub/affine.hh>
#include <dune/detailed/discretizations/evaluation/local/binary/elliptic.hh>
#include <dune/detailed/discretizations/discreteoperator/local/codim0/integral.hh>
#include <dune/detailed/discretizations/evaluation/local/unary/scale.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim0/integral.hh>
#include <dune/detailed/discretizations/discretefunctional/local/codim1/integral.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/matrix.hh>
#include <dune/detailed/discretizations/assembler/local/codim0/vector.hh>
#include <dune/detailed/discretizations/assembler/local/codim1/vector.hh>
#include <dune/detailed/discretizations/assembler/system.hh>


const std::string id = "elliptic.continuousgalerkin";

#ifndef POLORDER
  const int polOrder = 1;
#else
  const int polOrder = POLORDER;
#endif


/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(const std::string& filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "grid = " << "stuff.grid.provider.cube" << std::endl;
    file << "boundaryinfo = " << "stuff.grid.boundaryinfo.alldirichlet" << std::endl;
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = [12; 12; 12]" << std::endl;
    file << "[stuff.grid.boundaryinfo.idbased]" << std::endl;
    file << "dirichlet = [1; 2; 3]" << std::endl;
    file << "neumann = [4]" << std::endl;
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
    file << "expression = [0.1*x[0]; 0.0; 0.0]" << std::endl;
    file << "[neumann]" << std::endl;
    file << "order = 0"  << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 0.0; 0.0]" << std::endl;
    file << "[solver]" << std::endl;
    file << "type = bicgstab.diagonal"  << std::endl;
    file << "maxIter = 5000"  << std::endl;
    file << "precision = 1e-12"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()


int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".description";
    ensureParamFile(paramFilename);
    Dune::Stuff::Common::ExtendedParameterTree description(argc, argv, paramFilename);
    description.assertSub(id);

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                         Dune::Stuff::Common::LOG_CONSOLE);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();

    // timer
    Dune::Timer timer;

    info << "setting up grid:" << std::endl;
    typedef Dune::Stuff::GridProviderInterface<> GridProviderType;
    const GridProviderType* gridProvider
        = Dune::Stuff::GridProviders<>::create(description.get(id + ".grid", "stuff.grid.provider.cube"),
                                               description);
    typedef GridProviderType::GridType GridType;
    const Dune::shared_ptr< const GridType > grid = gridProvider->grid();
    typedef Dune::grid::Part::Leaf::Const< GridType > GridPartType;
    const GridPartType gridPart(*grid);
    typedef typename Dune::Stuff::Grid::BoundaryInfo::Interface< typename GridPartType::GridViewType > BoundaryInfoType;
    const Dune::shared_ptr< const BoundaryInfoType > boundaryInfo(
          Dune::Stuff::Grid::BoundaryInfo::create< typename GridPartType::GridViewType >(
              description.get(id + ".boundaryinfo", "stuff.grid.boundaryinfo.alldirichlet"),
                            description));

    info << "  took " << timer.elapsed() << " sec, has " << grid->size(0) << " entities" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(description.get(id + ".filename", id) + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    info << "initializing function space and data functions... " << std::flush;
    timer.reset();
    const int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const int DUNE_UNUSED(dimRange) = 1;
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    typedef Dune::FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    timer.reset();
    typedef Dune::Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange >
        ExpressionFunctionType;
    const Dune::shared_ptr< const ExpressionFunctionType >
        diffusion(ExpressionFunctionType::create(description.sub("diffusion")));
    const Dune::shared_ptr< const ExpressionFunctionType >
        force(ExpressionFunctionType::create(description.sub("force")));
    const Dune::shared_ptr< const ExpressionFunctionType >
        dirichlet(ExpressionFunctionType::create(description.sub("dirichlet")));
    const Dune::shared_ptr< const ExpressionFunctionType >
        neumann(ExpressionFunctionType::create(description.sub("neumann")));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "initializing discrete function spaces... " << std::flush;
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunctionSpace
        ::Continuous
        ::Lagrange< FunctionSpaceType, GridPartType, polOrder >
      LagrangeSpaceType;
    const LagrangeSpaceType lagrangeSpace(gridPart);
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunctionSpace
        ::Sub
        ::Linear
        ::Dirichlet< LagrangeSpaceType >
      TestSpaceType;
    const TestSpaceType testSpace(lagrangeSpace, boundaryInfo);
    typedef TestSpaceType AnsatzSpaceType;
    const AnsatzSpaceType ansatzSpace(lagrangeSpace, boundaryInfo);
    typedef typename Dune::Detailed::Discretizations::LA::Container::Factory::Eigen< RangeFieldType > ContainerFactory;
    typedef typename ContainerFactory::DenseVectorType VectorType;
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunction
        ::Default< LagrangeSpaceType, VectorType >
      DiscreteFunctionType;
    DiscreteFunctionType discreteDirichlet(lagrangeSpace, "dirichlet");
    Dune::Stuff::DiscreteFunction::Projection::Dirichlet::project(*boundaryInfo, *dirichlet, discreteDirichlet);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "initializing operator and functionals... " << std::flush;
    timer.reset();
    // * left hand side
    //   * elliptic operator
    typedef Dune::Detailed::Discretizations
        ::Evaluation
        ::Local
        ::Binary
        ::Elliptic< FunctionSpaceType, ExpressionFunctionType >
      EllipticEvaluationType;
    const EllipticEvaluationType ellipticEvaluation(diffusion, description.sub("diffusion").get("order", 0));
    typedef Dune::Detailed::Discretizations
        ::DiscreteOperator
        ::Local
        ::Codim0
        ::Integral< EllipticEvaluationType >
      EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(ellipticEvaluation);
    // * right hand side
    //   * L2 force functional
    typedef Dune::Detailed::Discretizations
        ::Evaluation
        ::Local
        ::Unary
        ::Scale< FunctionSpaceType, ExpressionFunctionType >
      ProductEvaluationType;
    const ProductEvaluationType forceEvaluation(force, description.sub("force").get("order", 0));
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunctional
        ::Local
        ::Codim0
        ::Integral< ProductEvaluationType >
      L2VolumeFunctionalType;
    const L2VolumeFunctionalType forceFunctional(forceEvaluation);
    //   * L2 neumann functional
    const ProductEvaluationType neumannEvaluation(neumann, description.sub("neumann").get("order", 0));
    typedef typename Dune::Detailed::Discretizations
        ::DiscreteFunctional
        ::Local
        ::Codim1
        ::Integral
        ::Boundary< ProductEvaluationType >
      L2BoundaryFunctionalType;
    const L2BoundaryFunctionalType neumannFunctional(neumannEvaluation);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "initializing matrix (of size " << testSpace.map().size() << "x" << ansatzSpace.map().size()
         << ") and vectors... " << std::flush;
    timer.reset();
    typedef ContainerFactory::RowMajorSparseMatrixType MatrixType;
    Dune::shared_ptr< MatrixType > systemMatrix = ContainerFactory::createRowMajorSparseMatrix(testSpace, ansatzSpace);
    Dune::shared_ptr< VectorType > forceVector = ContainerFactory::createDenseVector(testSpace);
    Dune::shared_ptr< VectorType > neumannVector = ContainerFactory::createDenseVector(testSpace);
    Dune::shared_ptr< VectorType > rhsVector = ContainerFactory::createDenseVector(testSpace);
    Dune::shared_ptr< VectorType > solutionVector = ContainerFactory::createDenseVector(testSpace);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "assembing system... " << std::flush;
    timer.reset();
    // * local matrix assembler
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Matrix< EllipticOperatorType >
        LocalMatrixAssemblerType;
    const Dune::shared_ptr< const LocalMatrixAssemblerType > localMatrixAssembler(
          new LocalMatrixAssemblerType(ellipticOperator));
    // * local vector assemblers
    //   * force vector
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim0::Vector< L2VolumeFunctionalType >
        LocalVolumeVectorAssemblerType;
    const Dune::shared_ptr< const LocalVolumeVectorAssemblerType > localforceVectorAssembler(
          new LocalVolumeVectorAssemblerType(forceFunctional));
    //   * neumann vector
    typedef Dune::Detailed::Discretizations::Assembler::Local::Codim1::Vector::Neumann< L2BoundaryFunctionalType,
                                                                                        BoundaryInfoType >
        LocalNeumannVectorAssemblerType;
    const Dune::shared_ptr< const LocalNeumannVectorAssemblerType > localNeumannVectorAssembler(
          new LocalNeumannVectorAssemblerType(neumannFunctional, boundaryInfo));
    typedef Dune::Detailed::Discretizations::Assembler::System< TestSpaceType, AnsatzSpaceType > SystemAssemblerType;
    // * system assembler
    SystemAssemblerType systemAssembler(testSpace, ansatzSpace);
    systemAssembler.addLocalMatrixAssembler(localMatrixAssembler, systemMatrix);
    systemAssembler.addLocalVectorAssembler(localforceVectorAssembler, forceVector);
    systemAssembler.addLocalVectorAssembler(localNeumannVectorAssembler, neumannVector);
    systemAssembler.assemble();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "applying constraints... " << std::flush;
    timer.reset();
    rhsVector->backend() = forceVector->backend()
        + neumannVector->backend()
        - systemMatrix->backend() * discreteDirichlet.vector()->backend();
    systemAssembler.applyConstraints(*systemMatrix, *rhsVector);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving linear system (of size " << systemMatrix->rows() << "x" << systemMatrix->cols() << ")" << std::endl;
    const std::string solverType = description.get("solver.type", "bicgstab");
    const unsigned int solverMaxIter = description.get("solver.maxIter", 5000);
    const double solverPrecision = description.get("solver.precision", 1e-12);
    info << "  using '" << solverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface< MatrixType, VectorType > SolverType;
    Dune::shared_ptr< SolverType > solver(Dune::Stuff::LA::Solver::create< MatrixType, VectorType >(solverType));
    const unsigned int failure = solver->apply(*systemMatrix,
                                               *rhsVector,
                                               *solutionVector,
                                               solverMaxIter,
                                               solverPrecision);
    if (failure)
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << solverType << "' reported a problem!");
    if (solutionVector->size() != ansatzSpace.map().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << solverType << "' produced a solution of wrong size (is "
                 << solutionVector->size() << ", should be " << ansatzSpace.map().size() << ")!");
    solutionVector->backend() += discreteDirichlet.vector()->backend();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const std::string solutionFilename = description.get(id + ".filename", id) + ".solution";
    const std::string solutionName = id + ".solution";
    info << "writing solution to '" << solutionFilename;
    if (dimDomain == 1)
      info << ".vtp";
    else
      info << ".vtu";
    info << "'... " << std::flush;
    timer.reset();
    typedef Dune::Detailed::Discretizations
        ::DiscreteFunction
        ::DefaultConst< LagrangeSpaceType, VectorType >
      ConstDiscreteFunctionType;
    const Dune::shared_ptr< const ConstDiscreteFunctionType > solution(new ConstDiscreteFunctionType(lagrangeSpace,
                                                                                                     solutionVector,
                                                                                                     solutionName));
    typedef Dune::VTKWriter< LagrangeSpaceType::GridViewType > VTKWriterType;
    VTKWriterType vtkWriter(lagrangeSpace.gridView());
    vtkWriter.addVertexData(solution);
    vtkWriter.write(solutionFilename);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
