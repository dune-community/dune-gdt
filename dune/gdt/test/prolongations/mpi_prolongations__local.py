# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017 - 2018)
#   Rene Milk       (2017 - 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

import spaces as sp
from grids import LevelGrids

dg = sp.DG(cache, base=LevelGrids)
fv = sp.FV(cache, base=LevelGrids, rdim=[1])
rt = sp.RT(cache, base=LevelGrids)
spaces = []
names = []
for sp, nm in zip(dg.spaces + fv.spaces + rt.spaces, dg.names + fv.names + rt.names):
    if not ('RaviartThomasSpace' in sp and 'Alu' in sp and 'LevelGridView' in sp):
        spaces.append(sp)
        names.append(nm)
spaces_with_names = zip(spaces, names)
