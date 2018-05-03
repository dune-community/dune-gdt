# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017 - 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

import spaces as sp

dg = sp.DG(cache)
fv = sp.FV(cache, rdim=[1])
rt = sp.RT(cache)
spaces = dg.spaces + fv.spaces + rt.spaces
names = dg.names + fv.names + rt.names
spaces_with_names = zip(spaces, names)
