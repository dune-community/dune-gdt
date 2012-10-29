#include "config.h"

// system
#include <vector>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/grid/part/leaf.hh>

// dune-fem
#include <dune/fem/misc/mpimanager.hh>

// dune-stuff
#include <dune/stuff/grid/provider/cube.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/discretefunctionspace/continuous/lagrange.hh>

const std::string id = "basefunctionset.continuous.lagrange";

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

template <class DiscreteFunctionSpaceType>
void inspect_basefunctions(const DiscreteFunctionSpaceType& space)
{
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  const GridPartType& gridPart                        = space.gridPart();
  const typename GridPartType::IndexSetType& indexSet = gridPart.indexSet();
  // walk the grid
  for (typename GridPartType::template Codim<0>::IteratorType entityIt = gridPart.template begin<0>();
       entityIt != gridPart.template end<0>();
       ++entityIt) {
    // entity
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
    const EntityType& entity                     = *entityIt;
    const typename EntityType::Geometry geometry = entity.geometry();
    typedef typename EntityType::Geometry::GlobalCoordinate CoordType;
    std::cout << "  entity " << indexSet.index(entity) << ", " << entity.type();
    // basefunctionset
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType BaseFunctionSetType;
    const BaseFunctionSetType baseFunctionSet = space.baseFunctionSet().local(entity);
    std::cout << ", " << baseFunctionSet.size() << " basis functions" << std::endl;
    // walk corners
    for (int c = 0; c < geometry.corners(); ++c) {
      const CoordType corner = geometry.corner(c);
      const CoordType local  = geometry.local(corner);
      std::cout << "    corner " << c << ", (" << corner << ")" << std::endl;
      typedef typename BaseFunctionSetType::RangeType RangeType;
      typedef typename BaseFunctionSetType::JacobianRangeType JacobianRangeType;
      std::vector<RangeType> evals(baseFunctionSet.size(), RangeType(0));
      std::vector<JacobianRangeType> jacobians(baseFunctionSet.size(), JacobianRangeType(0));
      baseFunctionSet.evaluate(local, evals);
      baseFunctionSet.jacobian(local, jacobians);
      // loop over basis functions
      for (unsigned int i = 0; i < baseFunctionSet.size(); ++i) {
        std::cout << "    - bf_" << i << ": evaluate(" << corner << ") = " << evals[i] << std::endl;
        std::cout << "    - bf_" << i << ": jacobian(" << corner << ") = " << jacobians[i] << std::endl;
      } // loop over basis functions
    } // walk corners
  } // walk the grid
} // void inspect_basefunctions(...)

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIManager::initialize(argc, argv);

    // timer
    Dune::Timer timer;

    // parameter
    Dune::ParameterTree paramTree;
    paramTree["level"] = "2";

    std::cout << "setting up grid:" << std::endl;
    timer.reset();
    typedef Dune::Stuff::Grid::Provider::UnitCube<> GridProviderType;
    const GridProviderType gridProvider(paramTree);
    typedef GridProviderType::GridType GridType;
    const GridType& grid = gridProvider.grid();
    typedef Dune::grid::Part::Leaf::Const<GridType> GridPartType;
    const GridPartType gridPart(grid);
    std::cout << "  took " << timer.elapsed() << " sec, has " << grid.size(0) << " entities" << std::endl;
    std::cout << "visualizing grid... " << std::flush;
    gridProvider.visualize(id + ".grid");
    std::cout << " done (took " << timer.elapsed() << " sek)" << std::endl;

    std::cout << "initializing spaces... " << std::flush;
    timer.reset();
    const int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const int DUNE_UNUSED(dimRange) = 1;
    typedef double DomainFieldType;
    typedef double RangeFieldType;
    typedef Dune::FunctionSpace<DomainFieldType, RangeFieldType, dimDomain, dimRange> FunctionSpaceType;
    typedef Dune::Detailed::Discretizations::DiscreteFunctionSpace::Continuous::
        Lagrange<FunctionSpaceType, GridPartType, polOrder> DiscreteFunctionSpaceType;
    const DiscreteFunctionSpaceType space(gridPart);
    std::cout << "done (took " << timer.elapsed() << " sec)" << std::endl;

    std::cout << "inspecting basefunctions:" << std::endl;
    inspect_basefunctions(space);

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
