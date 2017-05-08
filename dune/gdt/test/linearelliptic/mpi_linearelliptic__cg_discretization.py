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

from dune.xt.codegen import typeid_to_typedef_name, is_found


grids = ['Yasp2Grid']
try:
    if cache['dune-alugrid']:
        grids.extend(['AluSimplex2dGridType'])
        # 'AluConform2dGridType', disabled due to https://gitlab.dune-project.org/extensions/dune-alugrid/issues/17
except KeyError:
    pass
casenames = ['AO2013TestCase', 'ER2007TestCase', 'ESV2007TestCase', 'MixedBoundaryTestCase', 'Spe10Model1TestCase']
space_backends = ['fem', 'pdelab']
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]


if is_found(cache, 'dune-istl_DIR'):
    permutations = itertools.product(testcases, space_backends, ['istl_sparse'])
    permutations = [(t,s,l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations]
else:
    permutations = []
