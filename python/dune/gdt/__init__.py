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


def init_logger(max_info_level=-1,
                max_debug_level=-1,
                enable_warnings=True,
                enable_colors=True,
                info_color='blue',
                debug_color='darkgray',
                warning_color='red'):
    from ._gdt import init_logger as _init_logger
    initializers = [_init_logger]
    for module_name in ('xt.common', 'xt.grid', 'xt.functions', 'xt.la'):
        try:
            mm = import_module('dune.{}'.format(module_name))
            initializers.append(mm.init_logger)
        except ModuleNotFoundError:
            pass
    for initializer in initializers:
        initializer(max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color,
                        warning_color)


modules = ['spaces', 'assembler', 'discretefunction', 'projections']
#for space_backend in ('fem', 'pdelab'):
for space_backend in ('fem',):
    #for la_backend in ('common', 'eigen', 'istl'):
    for la_backend in ('istl',):
        modules.append('operators_elliptic_{}_{}'.format(space_backend, la_backend))
        #for grid in ('alberta', 'alu', 'ug', 'yasp'):
        for grid in ('alu', 'yasp'):
            modules.append('operators_elliptic_ipdg_{}_{}_{}'.format(grid, space_backend, la_backend))
            modules.append('functionals_l2_{}_{}_{}'.format(grid, space_backend, la_backend))
modules.append('bindings')

for module_name in modules:
    mod = import_module('.__{}'.format(module_name), 'dune.gdt')
    to_import = [name for name in mod.__dict__ if not name.startswith('_')]
    globals().update({name: mod.__dict__[name] for name in to_import})

