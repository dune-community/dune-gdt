# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2017)
#   Rene Milk       (2017 - 2018)
#   Tobias Leibner  (2017)


dune_pybindxi_add_module(__discretefunction EXCLUDE_FROM_ALL
                            discretefunction.bindings.cc)

dune_pybindxi_add_module(__projections EXCLUDE_FROM_ALL
                            projections.bindings.cc)

dune_pybindxi_add_module(__spaces EXCLUDE_FROM_ALL
                            spaces.bindings.cc)

dune_pybindxi_add_module(__assembler EXCLUDE_FROM_ALL
                            assembler.bindings.cc)

dune_pybindxi_add_module(__spaces_block EXCLUDE_FROM_ALL
                            playground/spaces/block.bindings.cc)

# the local operators
# * diffusive flux estimation
dune_pybindxi_add_module(__local_diffusive_flux_estimation_operator EXCLUDE_FROM_ALL
                            local/diffusive-flux-estimation-operator.bindings.cc)

# * elliptic ipdg
dune_pybindxi_add_module(__local_elliptic_ipdg_operators EXCLUDE_FROM_ALL
                            local/elliptic-ipdg-operators.bindings.cc)

# the operators
# * elliptic
dune_pybindxi_add_module(__operators_elliptic EXCLUDE_FROM_ALL
                            operators/elliptic.bindings.cc)
# * elliptic-ipdg
dune_pybindxi_add_module(__operators_elliptic_ipdg EXCLUDE_FROM_ALL
                            operators/elliptic-ipdg.bindings.cc)
# * flux reconstruction
dune_pybindxi_add_module(__operators_fluxreconstruction EXCLUDE_FROM_ALL
                            operators/fluxreconstruction.bindings.cc)

# * oswald interpolation
dune_pybindxi_add_module(__operators_oswaldinterpolation EXCLUDE_FROM_ALL
                            operators/oswaldinterpolation.bindings.cc)

# * L2 operator
dune_pybindxi_add_module(__operators_l2 EXCLUDE_FROM_ALL
                            operators/l2.bindings.cc)

# * weighted L2 operator
dune_pybindxi_add_module(__operators_weighted_l2 EXCLUDE_FROM_ALL
                            operators/weighted-l2.bindings.cc)

# * the operators from ESV2007
dune_pybindxi_add_module(__operators_ESV2007 EXCLUDE_FROM_ALL
                            playground/operators/ESV2007.bindings.cc)

# * the operators from OS2015
dune_pybindxi_add_module(__operators_OS2015 EXCLUDE_FROM_ALL
                            playground/operators/OS2015.bindings.cc)

# the functionals
# * the operators from RS2017
dune_pybindxi_add_module(__operators_RS2017 EXCLUDE_FROM_ALL
                            playground/operators/RS2017.bindings.cc)

dune_pybindxi_add_module(__functionals_elliptic_ipdg EXCLUDE_FROM_ALL
                            functionals/elliptic-ipdg.bindings.cc)

dune_pybindxi_add_module(__functionals_l2 EXCLUDE_FROM_ALL
                            functionals/l2.bindings.cc)

dune_pybindxi_add_module(__prolongations EXCLUDE_FROM_ALL
                            prolongations/prolongations.bindings.cc)
