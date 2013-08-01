// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

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
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/common/color.hh>

#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/la/containerfactory/eigen.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/swipdg-fluxes.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localfunctional/codim0.hh>
#include <dune/gdt/localfunctional/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>

const std::string id = "elliptic.discontinuousgalerkin-swip";

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif
dune_static_assert((polOrder > 0), "ERROR: polOrder hast to be positive!");

using namespace Dune::GDT;


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
    file << "grid = gridprovider.cube" << std::endl;
    file << "boundaryinfo = boundaryinfo.normalbased" << std::endl;
    file << "penaltyFactor = 10.0" << std::endl;
    file << "[gridprovider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = [12; 12; 12]" << std::endl;
    file << "[boundaryinfo.normalbased]" << std::endl;
    file << "default = dirichlet" << std::endl;
    file << "compare_tolerance = 1e-10" << std::endl;
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
} // void ensureSettings()


int main(int argc, char** argv)
{
  try {
    std::cout << Dune::Stuff::Common::colorString("WARNING:")
              << " there still is something wrong for nonzero dirichlet values!" << std::endl;

    // mpi
    Dune::Fem::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".settings";
    ensureSettings(paramFilename);
    Dune::Stuff::Common::ExtendedParameterTree settings(argc, argv, paramFilename);
    settings.assertSub(id);
    const double penaltyFactor = settings.get<double>(id + ".penaltyFactor");
    assert(penaltyFactor > 0.0 && "It better be!");

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO | Dune::Stuff::Common::LOG_CONSOLE);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
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

    info << "initializing space... " << std::flush;
    timer.reset();
    typedef DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GridPartType, polOrder, RangeFieldType, dimRange>
        SpaceType;
    const SpaceType space(gridPart);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // left hand side
    typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<ExpressionFunctionType>> EllipticOperatorType;
    const EllipticOperatorType ellipticOperator(*diffusion);
    typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDGFluxes::CouplingPrimal<ExpressionFunctionType>>
        CouplingOperatorType;
    const CouplingOperatorType ipdgCouplingOperator(*diffusion, penaltyFactor);
    typedef LocalOperator::
        Codim1BoundaryIntegral<LocalEvaluation::SWIPDGFluxes::BoundaryDirichletLHS<ExpressionFunctionType>>
            DirichletOperatorType;
    const DirichletOperatorType ipdgDirichletOperator(*diffusion, penaltyFactor);
    // right hand side
    typedef LocalFunctional::Codim0Integral<LocalEvaluation::Product<ExpressionFunctionType>> ForceFunctionalType;
    const ForceFunctionalType forceFunctional(*force);
    typedef LocalFunctional::Codim1Integral<LocalEvaluation::SWIPDGFluxes::BoundaryDirichletRHS<ExpressionFunctionType,
                                                                                                ExpressionFunctionType>>
        DirichletFunctionalType;
    const DirichletFunctionalType dirichletFunctional(*diffusion, *dirichlet, penaltyFactor);
    typedef LocalFunctional::Codim1Integral<LocalEvaluation::Product<ExpressionFunctionType>> NeumannFunctionalType;
    const NeumannFunctionalType neumannFunctional(*neumann);

    info << "initializing matrix (of size " << space.mapper().size() << "x" << space.mapper().size()
         << ") and vectors... " << std::flush;
    // timer.reset();
    typedef ContainerFactoryEigen<RangeFieldType> ContainerFactory;
    typedef ContainerFactory::RowMajorSparseMatrixType MatrixType;
    typedef ContainerFactory::DenseVectorType VectorType;
    std::shared_ptr<MatrixType> systemMatrix(ContainerFactory::createRowMajorSparseMatrix(space, space));
    std::shared_ptr<VectorType> rhsVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> solutionVector(ContainerFactory::createDenseVector(space));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "assembing system... " << std::flush;
    timer.reset();
    // * local matrix assembler
    typedef LocalAssembler::Codim0Matrix<EllipticOperatorType> LocalVolumeMatrixAssemblerType;
    const LocalVolumeMatrixAssemblerType diffusionMatrixAssembler(ellipticOperator);
    typedef LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> LocalCouplingMatrixAssemblerType;
    const LocalCouplingMatrixAssemblerType couplingMatrixAssembler(ipdgCouplingOperator);
    typedef LocalAssembler::Codim1BoundaryMatrix<DirichletOperatorType> LocalDirichletMatrixAssemblerType;
    const LocalDirichletMatrixAssemblerType dirichletMatrixAssembler(ipdgDirichletOperator);
    // * local vector assemblers
    typedef LocalAssembler::Codim0Vector<ForceFunctionalType> LocalForceVectorAssemblerType;
    const LocalForceVectorAssemblerType forceVectorAssembler(forceFunctional);
    typedef LocalAssembler::Codim1Vector<DirichletFunctionalType> LocalDirichletVectorAssembler;
    const LocalDirichletVectorAssembler dirichletVectorAssembler(dirichletFunctional);
    typedef LocalAssembler::Codim1Vector<NeumannFunctionalType> LocalNeumannVectorAssembler;
    const LocalNeumannVectorAssembler neumannVectorAssembler(neumannFunctional);
    // * system assembler
    typedef SystemAssembler<SpaceType, SpaceType> SystemAssemblerType;
    SystemAssemblerType systemAssembler(space);
    systemAssembler.addLocalAssembler(diffusionMatrixAssembler, *systemMatrix);
    systemAssembler.addLocalAssembler(
        couplingMatrixAssembler, SystemAssemblerType::AssembleOnInnerPrimally(), *systemMatrix);
    systemAssembler.addLocalAssembler(
        dirichletMatrixAssembler, SystemAssemblerType::AssembleOnDirichlet(*boundaryInfo), *systemMatrix);
    systemAssembler.addLocalAssembler(forceVectorAssembler, *rhsVector);
    systemAssembler.addLocalAssembler(
        dirichletVectorAssembler, SystemAssemblerType::AssembleOnDirichlet(*boundaryInfo), *rhsVector);
    systemAssembler.addLocalAssembler(
        neumannVectorAssembler, SystemAssemblerType::AssembleOnNeumann(*boundaryInfo), *rhsVector);
    systemAssembler.assemble();
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
    typedef DiscreteFunctionDefaultConst<SpaceType, VectorType> DiscreteFunctionType;
    const auto solution = std::make_shared<const DiscreteFunctionType>(space, solutionVector, solutionName);
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
