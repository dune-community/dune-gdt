import spaces as sp
from grids import LevelGrids

cg = sp.CG(cache, base=LevelGrids)
dg = sp.DG(cache, base=LevelGrids)
fv = sp.FV(cache, base=LevelGrids, rdim=[1])
rt = sp.RT(cache, base=LevelGrids)
spaces = cg.spaces + dg.spaces + fv.spaces + rt.spaces
names = cg.names + dg.names + fv.names + rt.names
spaces_with_names = zip(spaces, names)
