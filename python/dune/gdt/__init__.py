
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


from ._gdt import *
