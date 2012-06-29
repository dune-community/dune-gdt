
#include "config.h"

// system
#include <iostream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/parametertree.hh>

// dune-helper-tools
#include <dune/helper-tools/common/parametertree.hh>

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
    file << "visualize.grid = rb_grid_provider_cube_grid" << std::endl;
    file << "visualize.msGrid = rb_grid_provider_cube_msGrid" << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);
    // parameter
    const std::string filename = "continuous_galerkin.param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::HelperTools::Common::ParameterTree::init(argc, argv, filename);
    paramTree.report();

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
