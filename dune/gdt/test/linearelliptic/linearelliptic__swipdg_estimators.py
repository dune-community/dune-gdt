# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017 - 2018)

import itertools
from dune.xt.codegen import typeid_to_typedef_name, la_backends, is_found

grids = []
try:
    if cache['dune-alugrid']:
        grids.extend(['AluSimplex2dGridType', 'AluConform2dGridType'])
except KeyError:
    pass

# Only dirichlet zero and no neumann test cases are allowed here, so neither ESV2007TestCase nor MixedBoundaryTestCase!
casenames = ['AO2013TestCase', 'ESV2007TestCase']
try:
    cache['DXT_DISABLE_LARGE_TESTS']
except KeyError:
    casenames.append('Spe10Model1TestCase')
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

permutations = itertools.product(testcases, ('gdt',), la_backends(cache))
permutations = [(t, s, l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations]
