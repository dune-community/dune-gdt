# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2017)

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_helper_lib(_gdt_bindings_helper
                               assembler/system.bindings.cc
                               discretefunction/default.bindings.cc
                               #operators/elliptic-ipdg.bindings/alu_common_fem.cc
                               #operators/elliptic-ipdg.bindings/alu_common_gdt.cc
                               #operators/elliptic-ipdg.bindings/alu_common_pdelab.cc
                               #operators/elliptic-ipdg.bindings/alu_eigen_fem.cc
                               #operators/elliptic-ipdg.bindings/alu_eigen_gdt.cc
                               #operators/elliptic-ipdg.bindings/alu_eigen_pdelab.cc
                               operators/elliptic-ipdg.bindings/alu_istl_fem.cc
                               #operators/elliptic-ipdg.bindings/alu_istl_gdt.cc
                               #operators/elliptic-ipdg.bindings/alu_istl_pdelab.cc
                               #operators/elliptic-ipdg.bindings/yasp_common_fem.cc
                               #operators/elliptic-ipdg.bindings/yasp_common_gdt.cc
                               #operators/elliptic-ipdg.bindings/yasp_common_pdelab.cc
                               #operators/elliptic-ipdg.bindings/yasp_eigen_fem.cc
                               #operators/elliptic-ipdg.bindings/yasp_eigen_gdt.cc
                               #operators/elliptic-ipdg.bindings/yasp_eigen_pdelab.cc
                               #operators/elliptic-ipdg.bindings/yasp_istl_fem.cc
                               #operators/elliptic-ipdg.bindings/yasp_istl_gdt.cc
                               #operators/elliptic-ipdg.bindings/yasp_istl_pdelab.cc
                               #operators/elliptic.bindings/common_fem.cc
                               #operators/elliptic.bindings/common_gdt.cc
                               #operators/elliptic.bindings/common_pdelab.cc
                               #operators/elliptic.bindings/eigen_fem.cc
                               #operators/elliptic.bindings/eigen_gdt.cc
                               #operators/elliptic.bindings/eigen_pdelab.cc
                               operators/elliptic.bindings/istl_fem.cc
                               #operators/elliptic.bindings/istl_gdt.cc
                               #operators/elliptic.bindings/istl_pdelab.cc
                               projections/dirichlet.bindings.cc
                               spaces/cg.bindings.cc
                               spaces/dg.bindings.cc
                               #spaces/fv.bindings.cc
                               )
  set(pybindxi_modules)
    dune_pybindxi_add_module(_gdt bindings.cc)
    list(APPEND pybindxi_modules _gdt)
  foreach(_target assembler discretefunction projections spaces)
    dune_pybindxi_add_module(${_target} ${_target}.bindings.cc)
    list(APPEND pybindxi_modules ${_target})
  endforeach()
  foreach(_space_backend fem pdelab)
    foreach(_la_backend common eigen istl)
      dune_pybindxi_add_module(operators_elliptic_${_space_backend}_${_la_backend}
                               operators/elliptic.bindings/${_space_backend}_${_la_backend}.cc)
      list(APPEND pybindxi_modules operators_elliptic_${_space_backend}_${_la_backend})
      foreach(_grid alu alberta ug yasp)
        dune_pybindxi_add_module(operators_elliptic_ipdg_${_grid}_${_space_backend}_${_la_backend}
                                 operators/elliptic-ipdg.bindings/${_grid}_${_space_backend}_${_la_backend}.cc)
        list(APPEND pybindxi_modules operators_elliptic_ipdg_${_grid}_${_space_backend}_${_la_backend})
      endforeach()
    endforeach()
  endforeach()
  foreach(_target ${pybindxi_modules})
    target_link_dune_default_libraries(${_target})
    add_dune_alberta_flags(GRIDDIM 2 ${_target})
  endforeach()
  if(DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS ${pybindxi_modules})
  else()
    add_custom_target(bindings DEPENDS ${pybindxi_modules})
  endif()
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)

