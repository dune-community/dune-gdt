from dune.gdt.__spaces import *
from dune.gdt.__spaces_block import *

def make_block_dg_space(grid_provider):
    for factory in [globals()[s] for s in globals().keys() if s.startswith('make_block_dg')]:
        try:
            return factory(grid_provider)
        except:
            continue
    raise TypeError('no matching block dg space for {}'.format(grid_provider.__class__))

def make_dg_space(grid_provider):
    for factory in [globals()[s] for s in globals().keys() if s.startswith('make_dg')]:
        try:
            return factory(grid_provider)
        except:
            continue
    raise TypeError('no matching dg space for {}'.format(grid_provider.__class__))
