import itertools
from dune.xt.codegen import typeid_to_typedef_name, la_backends

grids = ['Yasp2Grid']
try:
    if cache['dune-alugrid']:
        grids.append('AluSimplex2dGridType')
        # 'AluConform2dGridType', diabled due to https://gitlab.dune-project.org/extensions/dune-alugrid/issues/17
except KeyError:
    pass
casenames = ['AO2013TestCase', 'ER2007TestCase', 'ESV2007TestCase', 'MixedBoundaryTestCase']
try:
    if bool(cache['DXT_ENABLE_LARGE_TESTS']):
        casenames.append('Spe10Model1TestCase')
except KeyError:
    pass

space_backends = ['fem', 'pdelab']

testcases = ['Dune::GDT::LinearElliptic::{}<{}>'.format(c, g) for c, g in itertools.product(casenames, grids)]

def filter(casename):
    if 'Spe10Model1TestCase' not in casename:
        return True
    return 'YaspGrid' in casename
permutations = itertools.product(testcases, space_backends, la_backends(cache))
permutations = [(t, s, l, typeid_to_typedef_name('{}_{}_{}'.format(t, s, l))) for t, s, l in permutations if filter(t)]
