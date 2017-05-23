# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)

import itertools
from dune.xt.codegen import typeid_to_typedef_name, la_backends

grids = []
try:
    if cache['dune-alugrid']:
        grids.extend(['AluConform2dGridType'])
except KeyError:
    pass

casenames = ['ESV2007DdSubdomainsTestCase',]
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

space_backends = []
for s in ('fem',):
    try:
        if cache['dune-{}'.format(s)]:
            space_backends.extend([s])
    except KeyError:
        pass

if len(space_backends) == 0 or len(la_backends(cache)) == 0:
    # prevent unusable iteration in template
    permutations = []
else:
    permutations = itertools.product(testcases, space_backends, ('istl_sparse', ))
    permutations = [(t, s, l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations]
