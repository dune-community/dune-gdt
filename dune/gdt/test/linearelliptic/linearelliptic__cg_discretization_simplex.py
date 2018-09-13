# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017 - 2018)
#   Rene Milk       (2017 - 2018)
#
# ~~~

import itertools

from dune.xt.codegen import typeid_to_typedef_name, la_backends, is_found


grid = []
# ///TODO re-add alugrid

casenames = ['AO2013TestCase', 'ER2007TestCase', 'ESV2007TestCase', 'MixedBoundaryTestCase',
             'Spe10Model1TestCase']
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

permutations = itertools.product(testcases, ('gdt',), la_backends(cache))
permutations = [(t,s,l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations]
