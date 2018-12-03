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

dune_pybindxi_add_module(__discretefunction EXCLUDE_FROM_ALL ${header} discretefunction/bindings.cc)

dune_pybindxi_add_module(__projections EXCLUDE_FROM_ALL ${header} projections/bindings.cc)

dune_pybindxi_add_module(__spaces EXCLUDE_FROM_ALL ${header} spaces/bindings.cc)

dune_pybindxi_add_module(__assembler
                         EXCLUDE_FROM_ALL
                         ${header}
                         assembler/system.cc
                         assembler/system/alu.cc
                         assembler/system/yasp.cc)

dune_pybindxi_add_module(__spaces_block EXCLUDE_FROM_ALL ${header} playground/spaces/block.cc)

# the local operators * diffusive flux estimation
dune_pybindxi_add_module(__local_diffusive_flux_estimation_operator
                         EXCLUDE_FROM_ALL
                         ${header}
                         local/diffusive-flux-estimation-operator.cc)

# * elliptic ipdg
dune_pybindxi_add_module(__local_elliptic_ipdg_operators EXCLUDE_FROM_ALL ${header} local/elliptic-ipdg-operators.cc)

# the operators * elliptic
dune_pybindxi_add_module(__operators_elliptic
                         EXCLUDE_FROM_ALL
                         ${header}
                         operators/elliptic/bindings.cc
                         operators/elliptic/istl.cc)
# * elliptic-ipdg
dune_pybindxi_add_module(__operators_elliptic_ipdg
                         EXCLUDE_FROM_ALL
                         ${header}
                         operators/elliptic-ipdg/bindings.cc
                         operators/elliptic-ipdg/alu_istl.cc
                         operators/elliptic-ipdg/yasp_istl.cc)
# * flux reconstruction
dune_pybindxi_add_module(__operators_fluxreconstruction EXCLUDE_FROM_ALL ${header} operators/fluxreconstruction.cc)

# * oswald interpolation
dune_pybindxi_add_module(__operators_oswaldinterpolation EXCLUDE_FROM_ALL ${header} operators/oswaldinterpolation.cc)

# * L2 operator
dune_pybindxi_add_module(__operators_l2 EXCLUDE_FROM_ALL ${header} operators/l2.cc)

# * weighted L2 operator
dune_pybindxi_add_module(__operators_weighted_l2 EXCLUDE_FROM_ALL ${header} operators/weighted-l2.cc)

# * the operators from ESV2007
dune_pybindxi_add_module(__operators_ESV2007 EXCLUDE_FROM_ALL ${header} playground/operators/ESV2007.cc)

# * the operators from OS2015
dune_pybindxi_add_module(__operators_OS2015 EXCLUDE_FROM_ALL ${header} playground/operators/OS2015.cc)

# the functionals * the operators from RS2017
dune_pybindxi_add_module(__operators_RS2017 EXCLUDE_FROM_ALL ${header} playground/operators/RS2017.cc)

dune_pybindxi_add_module(__functionals_elliptic_ipdg
                         EXCLUDE_FROM_ALL
                         ${header}
                         functionals/elliptic-ipdg/bindings.cc
                         functionals/elliptic-ipdg/alu_istl.cc
                         functionals/elliptic-ipdg/yasp_istl.cc)

dune_pybindxi_add_module(__functionals_l2
                         EXCLUDE_FROM_ALL
                         ${header}
                         functionals/l2/bindings.cc
                         functionals/l2/alu_istl.cc
                         functionals/l2/yasp_istl.cc)
