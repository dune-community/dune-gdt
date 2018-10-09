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

all = sp.all_spaces_with_names_and_grids(cache, base=LevelGrids, rdim=[1])
spaces_with_names = [(space, nm, grid) for space, nm, grid in all if 'AluConform2dLevelGrid' not in space]
