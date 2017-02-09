import itertools

from dune.xt.codegen import typeid_to_typedef_name


grids = ['Yasp2Grid']
try:
    if cache['dune-alugrid']:
        grids.extend(['AluConform2dGridType', 'AluSimplex2dGridType'])
except KeyError:
    pass
casenames = ['AO2013TestCase', 'ER2007TestCase', 'ESV2007TestCase', 'MixedBoundaryTestCase', 'Spe10Model1TestCase']
space_backends = ['fem', 'pdelab']
la_backends = ['eigen_sparse', 'istl_sparse']
testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

permutations = itertools.product(testcases, space_backends, la_backends)
permutations = [(t,s,l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations]