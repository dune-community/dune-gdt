# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

if(NOT CMAKE_WITH_AUTOTOOLS)

  set(BASEDIR ../..)
  set(CURDIR doc/doxygen)
  # Create Doxyfile.in and Doxyfile, and doxygen documentation
  add_doxygen_target()

endif(NOT CMAKE_WITH_AUTOTOOLS)
