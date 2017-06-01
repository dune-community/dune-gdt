# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017)

import itertools
from dune.xt.codegen import typeid_to_typedef_name, la_backends

grids = ['Yasp2Grid']
try:
    if cache['dune-alugrid']:
        grids.extend(['AluSimplex2dGridType', 'AluConform2dGridType'])
except KeyError:
    pass

casenames = ['AO2013TestCase', 'ER2007TestCase', 'ESV2007TestCase', 'MixedBoundaryTestCase']
try:
    cache['DXT_DISABLE_LARGE_TESTS']
except KeyError:
    casenames.append('Spe10Model1TestCase')
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

space_backends = []
for s in ('fem',):
    try:
        if cache['dune-{}'.format(s)]:
            space_backends.extend([s])
    except KeyError:
        pass

if len(space_backends) == 0:
    # prevent unusable iteration in template
    permutations = []
else:
    permutations = itertools.product(testcases, space_backends, la_backends(cache))

def filter(t, s):
    # pdelab has no DG impl for simplicial grids
    if s == 'pdelab':
        return 'AluSimplex' not in t
    return True

permutations = [(t, s, l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations if filter(t, s)]
