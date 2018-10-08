# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017 - 2018)
#
# ~~~

import spaces as sp
from grids import LevelGrids

fv = sp.FV(cache, base=LevelGrids, rdim=[1])
spaces_with_names = []
for sp, nm, grid in zip(fv.spaces ,fv.names, fv.grids ):
    if not ('AluConform2dLevelGrid' in sp):
        spaces_with_names.append((sp, nm, grid))
