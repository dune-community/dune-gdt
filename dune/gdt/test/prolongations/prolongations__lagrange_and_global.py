import spaces as sp
from grids import LevelGrids

cg = sp.CG(cache, base=LevelGrids)
spaces = cg.spaces
names = cg.names
spaces_with_names = zip(spaces, names)
