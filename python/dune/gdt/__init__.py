# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017 - 2018)
#   Rene Milk       (2016, 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

from importlib import import_module

import numpy as np

import dune.xt.common # inits MPI via mpi4py
import dune.xt.la
import dune.xt.grid
import dune.xt.functions

# The ordering of these module imports matter, do not change unless you know what you are doing!
_gdt_modules = ['spaces',
                'spaces_block',
                'local_diffusive_flux_estimation_operator',
                'local_elliptic_ipdg_operators',
                'assembler',
                'discretefunction',
                'projections',
                'prolongations',
                'functionals_elliptic_ipdg',
                'functionals_l2',
                'operators_elliptic',
                'operators_elliptic_ipdg',
                'operators_fluxreconstruction',
                'operators_oswaldinterpolation',
                'operators_ESV2007',
                'operators_OS2015',
                'operators_RS2017',
                'operators_l2',
                'operators_weighted_l2']

for module_name in _gdt_modules:
    import_module('.__{}'.format(module_name), 'dune.gdt')
