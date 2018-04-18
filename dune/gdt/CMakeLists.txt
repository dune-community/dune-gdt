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
  # some top-level modules
  foreach(_target assembler discretefunction projections prolongations spaces usercode)
    dune_pybindxi_add_module(__${_target} ${_target}.bindings.cc)
    list(APPEND pybindxi_modules __${_target})
  endforeach()
  dune_pybindxi_add_module(__spaces_block playground/spaces/block.bindings.cc)
  list(APPEND pybindxi_modules __spaces_block)
  # the local operators
  # * diffusive flux estimation
  dune_pybindxi_add_module(__local_diffusive_flux_estimation_operator local/diffusive-flux-estimation-operator.bindings.cc)
  list(APPEND pybindxi_modules __local_diffusive_flux_estimation_operator)
  # * elliptic ipdg
  dune_pybindxi_add_module(__local_elliptic_ipdg_operators local/elliptic-ipdg-operators.bindings.cc)
  list(APPEND pybindxi_modules __local_elliptic_ipdg_operators)
  # the operators
  # * elliptic
  dune_pybindxi_add_module(__operators_elliptic operators/elliptic.bindings.cc)
  list(APPEND pybindxi_modules __operators_elliptic)
  # * elliptic-ipdg
  dune_pybindxi_add_module(__operators_elliptic_ipdg operators/elliptic-ipdg.bindings.cc)
  list(APPEND pybindxi_modules __operators_elliptic_ipdg)
  # * flux reconstruction
  dune_pybindxi_add_module(__operators_fluxreconstruction operators/fluxreconstruction.bindings.cc)
  list(APPEND pybindxi_modules __operators_fluxreconstruction)
  # * oswald interpolation
  dune_pybindxi_add_module(__operators_oswaldinterpolation operators/oswaldinterpolation.bindings.cc)
  list(APPEND pybindxi_modules __operators_oswaldinterpolation)
  # * L2 operator
  dune_pybindxi_add_module(__operators_l2 operators/l2.bindings.cc)
  list(APPEND pybindxi_modules __operators_l2)
  # * weighted L2 operator
  dune_pybindxi_add_module(__operators_weighted_l2 operators/weighted-l2.bindings.cc)
  list(APPEND pybindxi_modules __operators_weighted_l2)
  # * the operators from ESV2007
  dune_pybindxi_add_module(__operators_ESV2007 playground/operators/ESV2007.bindings.cc)
  list(APPEND pybindxi_modules __operators_ESV2007)
  # * the operators from OS2015
  dune_pybindxi_add_module(__operators_OS2015 playground/operators/OS2015.bindings.cc)
  list(APPEND pybindxi_modules __operators_OS2015)
  # the functionals
  # * the operators from RS2017
  dune_pybindxi_add_module(__operators_RS2017 playground/operators/RS2017.bindings.cc)
  list(APPEND pybindxi_modules __operators_RS2017)
  dune_pybindxi_add_helper_lib(__functionals_elliptic_ipdg_helper
                               functionals/elliptic-ipdg.bindings/alu_fem_istl.cc
                               functionals/elliptic-ipdg.bindings/yasp_fem_istl.cc)
  dune_pybindxi_add_module(__functionals_elliptic_ipdg functionals/elliptic-ipdg.bindings.cc)
  target_link_libraries(__functionals_elliptic_ipdg __functionals_elliptic_ipdg_helper)
  list(APPEND pybindxi_modules __functionals_elliptic_ipdg)
  dune_pybindxi_add_module(__functionals_l2 functionals/l2.bindings.cc)
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

