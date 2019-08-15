# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2018)
# ~~~

from dune.xt import guarded_import

from ._discretefunction_1d import *
from ._discretefunction_2d import *
from ._discretefunction_3d import *


def make_discrete_function(space, vector, name="dune.gdt.discretefunction"):
    for factory in [globals()[s] for s in globals().keys() if s.startswith('make_discrete_function_')]:
        try:
            return factory(space, vector, name)
        except:
            continue
    raise TypeError('no matching factory for space \'{}\''.format(space))
