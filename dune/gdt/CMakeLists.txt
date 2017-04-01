# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2017)

if(dune-pybindxi_FOUND)
  set(pybindxi_modules)
    dune_pybindxi_add_module(__bindings bindings.cc)
    list(APPEND pybindxi_modules __bindings)
  foreach(_target assembler discretefunction projections spaces)
    dune_pybindxi_add_module(__${_target} ${_target}.bindings.cc)
    list(APPEND pybindxi_modules __${_target})
  endforeach()
  #foreach(_space_backend fem pdelab)
  foreach(_space_backend fem)
    #foreach(_la_backend common eigen istl)
    foreach(_la_backend istl)
      dune_pybindxi_add_module(__operators_elliptic_${_space_backend}_${_la_backend}
                               operators/elliptic.bindings/${_space_backend}_${_la_backend}.cc)
      list(APPEND pybindxi_modules __operators_elliptic_${_space_backend}_${_la_backend})
      #foreach(_grid alu alberta ug yasp)
      foreach(_grid alu yasp)
        dune_pybindxi_add_module(__operators_elliptic_ipdg_${_grid}_${_space_backend}_${_la_backend}
                                 operators/elliptic-ipdg.bindings/${_grid}_${_space_backend}_${_la_backend}.cc)
        list(APPEND pybindxi_modules __operators_elliptic_ipdg_${_grid}_${_space_backend}_${_la_backend})
      endforeach()
    endforeach()
  endforeach()
  dune_pybindxi_add_helper_lib(__functionals_l2_helper
                               functionals/l2.bindings/alu_fem_istl.cc
                               functionals/l2.bindings/yasp_fem_istl.cc)
  dune_pybindxi_add_module(__functionals_l2 functionals/l2.bindings.cc)
  target_link_libraries(__functionals_l2 __functionals_l2_helper)
  list(APPEND pybindxi_modules __functionals_l2)
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

