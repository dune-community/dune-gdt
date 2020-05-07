# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2018)
#   René Fritze     (2018)
#   Tobias Leibner  (2018)
# ~~~

file(GLOB_RECURSE header "*.hh")

dune_pybindxi_add_module(_discretefunction EXCLUDE_FROM_ALL ${header} discretefunction/discretefunction.cc)
dune_pybindxi_add_module(_local_integrands_binary_element_interface EXCLUDE_FROM_ALL ${header} local/integrands/binary_element_interface.cc)
dune_pybindxi_add_module(_local_integrands_laplace EXCLUDE_FROM_ALL ${header} local/integrands/laplace.cc)
dune_pybindxi_add_module(_operators_interfaces_common EXCLUDE_FROM_ALL ${header} operators/interfaces_common.cc)
dune_pybindxi_add_module(_operators_interfaces_eigen EXCLUDE_FROM_ALL ${header} operators/interfaces_eigen.cc)
dune_pybindxi_add_module(_operators_interfaces_istl EXCLUDE_FROM_ALL ${header} operators/interfaces_istl.cc)
dune_pybindxi_add_module(_spaces_h1_continuous_lagrange EXCLUDE_FROM_ALL ${header} spaces/h1/continuous_lagrange.cc)
dune_pybindxi_add_module(_spaces_interface EXCLUDE_FROM_ALL ${header} spaces/interface.cc)
# dune_pybindxi_add_module(usercode EXCLUDE_FROM_ALL ${header} usercode.cc)
# dune_pybindxi_add_module(gamm_2019_talk_on_conservative_rb
#                          EXCLUDE_FROM_ALL
#                          ${header}
#                          gamm-2019-talk-on-conservative-rb.cc)

