#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#define HAVE_DUNE_DETAILED_DISCRETIZATIONS 1

#include <iostream>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
#include <dune/common/dynvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>

#include <dune/detailed/discretizations/space/continuouslagrange/fem.hh>
#include <dune/detailed/discretizations/la/containerfactory/eigen.hh>
#include <dune/detailed/discretizations/localevaluation/elliptic.hh>
#include <dune/detailed/discretizations/localevaluation/product.hh>
#include <dune/detailed/discretizations/localoperator/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim1.hh>
#include <dune/detailed/discretizations/assembler/local/codim0.hh>
#include <dune/detailed/discretizations/assembler/local/codim1.hh>
#include <dune/detailed/discretizations/space/constraints.hh>
#include <dune/detailed/discretizations/assembler/system.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>


const std::string id = "elliptic.continuousgalerkin";

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif
dune_static_assert((polOrder > 0), "ERROR: polOrder hast to be positive!");

using namespace Dune::Detailed::Discretizations;


/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureSettings(const std::string& filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "grid = "
         << "gridprovider.cube" << std::endl;
    file << "boundaryinfo = "
         << "boundaryinfo.normalbased" << std::endl;
    file << "[gridprovider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = [12; 12; 12]" << std::endl;
    file << "[boundaryinfo.normalbased]" << std::endl;
    file << "default = dirichlet" << std::endl;
    file << "neumann = [0.0; 1.0]" << std::endl;
    file << "[diffusion]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[force]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[dirichlet]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [0.1*x[0]; 0.0; 0.0]" << std::endl;
    file << "[neumann]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [0.1; 0.0; 0.0]" << std::endl;
    file << "[solver]" << std::endl;
    file << "type = bicgstab.ilut" << std::endl;
    file << "maxIter = 5000" << std::endl;
    file << "precision = 1.0e-8" << std::endl;
    file << "preconditioner.dropTol = 1e-4" << std::endl;
    file << "preconditioner.fillFactor = 10" << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureSettings(...)


int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::Fem::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".settings";
    ensureSettings(paramFilename);
    Dune::Stuff::Common::ExtendedParameterTree settings(argc, argv, paramFilename);
    settings.assertSub(id);

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO | Dune::Stuff::Common::LOG_CONSOLE);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();

    // timer
    Dune::Timer timer;

    info << "setting up grid:" << std::endl;
    typedef Dune::Stuff::GridProviderInterface<> GridProviderType;
    const GridProviderType* gridProvider =
        Dune::Stuff::GridProviders<>::create(settings.get(id + ".grid", "gridprovider.cube"), settings);
    typedef GridProviderType::GridType GridType;
    const std::shared_ptr<const GridType> grid = gridProvider->grid();
    typedef Dune::grid::Part::Leaf::Const<GridType> GridPartType;
    typedef typename GridPartType::GridViewType GridViewType;
    const GridPartType gridPart(*grid);
    typedef typename Dune::Stuff::GridboundaryInterface<typename GridPartType::GridViewType> BoundaryInfoType;
    const std::shared_ptr<const BoundaryInfoType> boundaryInfo(
        Dune::Stuff::Gridboundaries<typename GridPartType::GridViewType>::create(
            settings.get(id + ".boundaryinfo", "boundaryinfo.alldirichlet"), settings));
    info << "  took " << timer.elapsed() << " sec, has " << grid->size(0) << " entities" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(settings.get(id + ".boundaryinfo", "boundaryinfo.alldirichlet"),
                            settings,
                            settings.get(id + ".filename", id) + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef typename GridType::ctype DomainFieldType;
    typedef double RangeFieldType;

    typedef Dune::Stuff::FunctionExpression<DomainFieldType, dimDomain, RangeFieldType, dimRange>
        ExpressionFunctionType;
    const std::shared_ptr<const ExpressionFunctionType> diffusion(
        ExpressionFunctionType::create(settings.sub("diffusion")));
    const std::shared_ptr<const ExpressionFunctionType> force(ExpressionFunctionType::create(settings.sub("force")));
    const std::shared_ptr<const ExpressionFunctionType> dirichlet(
        ExpressionFunctionType::create(settings.sub("dirichlet")));
    const std::shared_ptr<const ExpressionFunctionType> neumann(
        ExpressionFunctionType::create(settings.sub("neumann")));

    info << "initializing discrete function spaces... " << std::flush;
    timer.reset();
    typedef ContinuousLagrangeSpace::FemWrapper<GridPartType, polOrder, RangeFieldType, dimRange> SpaceType;
    const SpaceType space(gridPart);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // left hand side
    // * elliptic diffusion operator
    typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<ExpressionFunctionType>> EllipticOperatorType;
    const EllipticOperatorType diffusionOperator(*diffusion);
    // * right hand side
    //   * L2 force functional
    typedef LocalFunctional::Codim0Integral<LocalEvaluation::Product<ExpressionFunctionType>> L2VolumeFunctionalType;
    const L2VolumeFunctionalType forceFunctional(*force);
    //   * L2 neumann functional
    typedef LocalFunctional::Codim1Integral<LocalEvaluation::Product<ExpressionFunctionType>> L2FaceFunctionalType;
    const L2FaceFunctionalType neumannFunctional(*neumann);

    info << "initializing matrix (of size " << space.mapper().size() << "x" << space.mapper().size()
         << ") and vectors... " << std::flush;
    timer.reset();
    typedef ContainerFactoryEigen<RangeFieldType> ContainerFactory;
    typedef ContainerFactory::RowMajorSparseMatrixType MatrixType;
    typedef ContainerFactory::DenseVectorType VectorType;
    std::shared_ptr<MatrixType> systemMatrix(ContainerFactory::createRowMajorSparseMatrix(space, space));
    std::shared_ptr<VectorType> forceVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> dirichletVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> neumannVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> rhsVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> solutionVector(ContainerFactory::createDenseVector(space));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "assembing system... " << std::flush;
    timer.reset();
    // * dirichlet boundary values
    typedef DiscreteFunctionDefault<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType dirichletProjection(space, dirichletVector, "dirichlet");
    Dune::Stuff::DiscreteFunction::project(*boundaryInfo, *dirichlet, dirichletProjection);

    // * local matrix assembler
    typedef LocalAssembler::Codim0Matrix<EllipticOperatorType> LocalEllipticOperatorMatrixAssemblerType;
    const LocalEllipticOperatorMatrixAssemblerType diffusionMatrixAssembler(diffusionOperator);
    // * local vector assemblers
    //   * force vector
    typedef LocalAssembler::Codim0Vector<L2VolumeFunctionalType> LocalL2VolumeFunctionalVectorAssemblerType;
    const LocalL2VolumeFunctionalVectorAssemblerType forceVectorAssembler(forceFunctional);
    //   * neumann vector
    typedef LocalAssembler::Codim1Vector<L2FaceFunctionalType> LocalL2FaceFunctionalVectorAssemblerType;
    const LocalL2FaceFunctionalVectorAssemblerType neumannVectorAssembler(neumannFunctional);
    // * system assembler
    typedef SystemAssembler<SpaceType, SpaceType> SystemAssemblerType;
    SystemAssemblerType systemAssembler(space);
    systemAssembler.addLocalAssembler(diffusionMatrixAssembler, *systemMatrix);
    systemAssembler.addLocalAssembler(forceVectorAssembler, *forceVector);
    systemAssembler.addLocalAssembler(
        neumannVectorAssembler, SystemAssemblerType::AssembleOnNeumann(*boundaryInfo), *neumannVector);
    systemAssembler.assemble();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "applying constraints... " << std::flush;
    timer.reset();
    Constraints::Dirichlet<GridViewType, RangeFieldType> dirichletConstraints(
        *boundaryInfo, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
    rhsVector->backend() =
        forceVector->backend() + neumannVector->backend() - systemMatrix->backend() * dirichletVector->backend();
    systemAssembler.addLocalConstraints(dirichletConstraints, *systemMatrix);
    systemAssembler.addLocalConstraints(dirichletConstraints, *rhsVector);
    systemAssembler.applyConstraints();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving linear system (of size " << systemMatrix->rows() << "x" << systemMatrix->cols() << ")"
         << std::endl;
    const Dune::Stuff::Common::ExtendedParameterTree linearSolverSettings = settings.sub("solver");
    const std::string solverType                                          = linearSolverSettings.get("type", "bicgstab.ilut");
    info << "  using '" << solverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::SolverInterface<MatrixType, VectorType> SolverType;
    std::shared_ptr<SolverType> solver(Dune::Stuff::LA::createSolver<MatrixType, VectorType>(solverType));
    const size_t failure = solver->apply(*systemMatrix, *rhsVector, *solutionVector, linearSolverSettings);
    if (failure)
      DUNE_THROW(Dune::MathError, "\nERROR: linear solver '" << solverType << "' reported a problem!");
    if (solutionVector->size() != space.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << solverType << "' produced a solution of wrong size (is "
                                            << solutionVector->size()
                                            << ", should be "
                                            << space.mapper().size()
                                            << ")!");
    solutionVector->backend() += dirichletVector->backend();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const std::string solutionFilename = settings.get(id + ".filename", id) + ".solution";
    const std::string solutionName     = id + ".solution";
    info << "writing solution to '" << solutionFilename;
    if (dimDomain == 1)
      info << ".vtp";
    else
      info << ".vtu";
    info << "'... " << std::flush;
    timer.reset();
    typedef DiscreteFunctionDefaultConst<SpaceType, VectorType> ConstDiscreteFunctionType;
    const std::shared_ptr<const ConstDiscreteFunctionType> solution(
        new ConstDiscreteFunctionType(space, solutionVector, solutionName));
    typedef Dune::VTKWriter<GridViewType> VTKWriterType;
    VTKWriterType vtkWriter(gridPart.gridView());
    vtkWriter.addVertexData(solution);
    vtkWriter.write(solutionFilename);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

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
