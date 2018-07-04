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

cg = sp.CG(cache, base=LevelGrids)
dg = sp.DG(cache, base=LevelGrids)
spaces = []
names = []
for sp, nm in zip(cg.spaces + dg.spaces, cg.names + dg.names):
    if not ('AluConform2dLevelGrid' in sp):
        spaces.append(sp)
        names.append(nm)
spaces_with_names = zip(spaces, names)
