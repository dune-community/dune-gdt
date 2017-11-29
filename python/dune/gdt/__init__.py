# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2016)

from importlib import import_module

import numpy as np

from dune.xt.common import DEBUG # inits MPI via mpi4py
import dune.xt.la
import dune.xt.grid
import dune.xt.functions

_init_logger_methods = list()
_test_logger_methods = list()
_init_mpi_methods = list()
_other_modules = ('xt.common', 'xt.grid', 'xt.functions', 'xt.la')

# The ordering of these module imports matter, do not change unless you know what you are doing!
_gdt_modules = ['spaces',
                'spaces_block',
                'local_diffusive_flux_estimation_operator',
                'local_elliptic_ipdg_operators',
                'assembler',
                'discretefunction',
                'projections',
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
    mod = import_module('.__{}'.format(module_name), 'dune.gdt')
    to_import = [name for name in mod.__dict__ if not name.startswith('_')]
    globals().update({name: mod.__dict__[name] for name in to_import})
    _init_logger_methods.append(mod.__dict__['_init_logger'])
    _test_logger_methods.append(mod.__dict__['_test_logger'])
    _init_mpi_methods.append(mod.__dict__['_init_mpi'])

del _gdt_modules


def init_logger(max_info_level=999,
                max_debug_level=999,
                enable_warnings=True,
                enable_colors=True,
                info_color='blue',
                debug_color='darkgray',
                warning_color='red'):
    init_logger_methods = _init_logger_methods.copy()
    for module_name in _other_modules:
        try:
            mm = import_module('dune.{}'.format(module_name))
            for init_logger_method in mm._init_logger_methods:
                init_logger_methods.append(init_logger_method)
        except ModuleNotFoundError:
            pass
    for init_logger_method in init_logger_methods:
        init_logger_method(max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color,
                           warning_color)

def test_logger(info=True, debug=True, warning=True):
    test_logger_methods = _test_logger_methods.copy()
    for module_name in _other_modules:
        try:
            mm = import_module('dune.{}'.format(module_name))
            for test_logger_method in mm._test_logger_methods:
                test_logger_methods.append(test_logger_method)
        except ModuleNotFoundError:
            pass
    for test_logger_method in test_logger_methods:
        test_logger_method(info, debug, warning)

def init_mpi(args=list()):
    if DEBUG:
        init_mpi_methods = [_init_mpi_methods[0],]
    else:
        init_mpi_methods = _init_mpi_methods.copy()
        for module_name in _other_modules:
            try:
                mm = import_module('dune.{}'.format(module_name))
                for init_mpi_method in mm._init_mpi_methods:
                    init_mpi_methods.append(init_mpi_method)
            except ModuleNotFoundError:
                pass
    for init_mpi_method in init_mpi_methods:
        init_mpi_method(args)


HAVE_DUNE_FEM = np.any(['FemP1Space' in var for var in globals().keys()])
HAVE_DUNE_PDELAB = np.any(['PdelabP1Space' in var for var in globals().keys()])

